#include "RigidNURBS_1_RigidNURBS_1.h"


#include "RigidNURBS_1.h"
#include "SSContactData.h"
#include "NURBSSurface.h"
#include "Node.h"
#include "Dynamic.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

RigidNURBS_1_RigidNURBS_1::RigidNURBS_1_RigidNURBS_1()
{
	DefaultValues();
	derivative_order = 3;
	dataA = new Matrix*[derivative_order + 1];
	for (int i = 0; i < derivative_order + 1; i++)
	{
		dataA[i] = new Matrix[derivative_order + 1];
		for (int j = 0; j < derivative_order + 1; j++)
			dataA[i][j] = Matrix(3);
	}
	dataB = new Matrix*[derivative_order + 1];
	for (int i = 0; i < derivative_order + 1; i++)
	{
		dataB[i] = new Matrix[derivative_order + 1];
		for (int j = 0; j < derivative_order + 1; j++)
			dataB[i][j] = Matrix(3);
	}

	//Last calls of derivatives - improves efficiency
	last_cp = Matrix(4);
	last_ci = Matrix(4);
	first_ci = true;
	first_cp = true;
}


RigidNURBS_1_RigidNURBS_1::~RigidNURBS_1_RigidNURBS_1()
{
	Free();

	for (int i = 0; i < derivative_order + 1; i++)
		delete[] dataA[i];
	delete[]dataA;
	for (int i = 0; i < derivative_order + 1; i++)
		delete[] dataB[i];
	delete[]dataB;
}

//Chute inicial para coordenadas convectivas do par de superfícies
void RigidNURBS_1_RigidNURBS_1::InitialGuess(SSContactData* c_data)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	RigidNURBS_1* surfA;		//Ponteiro para a superfície A
	RigidNURBS_1* surfB;		//Ponteiro para a superfície B
	surfA = static_cast<RigidNURBS_1*>(db.surfaces[surf1_ID - 1]);
	surfB = static_cast<RigidNURBS_1*>(db.surfaces[surf2_ID - 1]);

	c_data->convective[0][0] = (surfA->surf->U_knot_vector[surfA->surf->U_order + surfA->surf->U_dim] + surfA->surf->U_knot_vector[surfA->surf->U_order]) / 2;
	c_data->convective[0][1] = (surfA->surf->V_knot_vector[surfA->surf->V_order + surfA->surf->V_dim] + surfA->surf->V_knot_vector[surfA->surf->V_order]) / 2;
	c_data->convective[0][2] = (surfB->surf->U_knot_vector[surfB->surf->U_order + surfB->surf->U_dim] + surfB->surf->U_knot_vector[surfB->surf->U_order]) / 2;
	c_data->convective[0][3] = (surfB->surf->V_knot_vector[surfB->surf->V_order + surfB->surf->V_dim] + surfB->surf->V_knot_vector[surfB->surf->V_order]) / 2;

	for (int ip = 1; ip < c_data->n_solutions; ip++)
	{
		//Preenchendo as coordenadas convectivas:
		c_data->convective[ip][0] = c_data->convective[0][0];
		c_data->convective[ip][1] = c_data->convective[0][1];
		c_data->convective[ip][2] = c_data->convective[0][2];
		c_data->convective[ip][3] = c_data->convective[0][3];
	}

}

//Evaluates the necessary derivatives to be used as input to AceGen routines
void RigidNURBS_1_RigidNURBS_1::EvaluateNURBSDerivatives_p(Matrix& mc)
{
	if (mc == last_cp && first_cp == false)
		return;
	else
	{
		last_cp = mc;
		first_cp = false;
	}
		
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	RigidNURBS_1* surfA;		//Ponteiro para a superfície A
	RigidNURBS_1* surfB;		//Ponteiro para a superfície B
	surfA = static_cast<RigidNURBS_1*>(db.surfaces[surf1_ID - 1]);
	surfB = static_cast<RigidNURBS_1*>(db.surfaces[surf2_ID - 1]);
	//Evaluates data and store in dataA and dataB
	surfA->surf->NURBSDerivatives(mc(0, 0), mc(1, 0), dataA, derivative_order);
	surfB->surf->NURBSDerivatives(mc(2, 0), mc(3, 0), dataB, derivative_order);
	//Copies data from dataA and dataB to AceGen input data pointers
	cAp[0] = mc(0, 0);
	cAp[1] = mc(1, 0);
	cBp[0] = mc(2, 0);
	cBp[1] = mc(3, 0);
	for (int i = 0; i < 3; i++)
	{
		GAp[i] = dataA[0][0](i, 0);
		dGAp[i][0] = dataA[1][0](i, 0);
		dGAp[i][1] = dataA[0][1](i, 0);
		ddGAp[i][0][0] = dataA[2][0](i, 0);
		ddGAp[i][0][1] = dataA[1][1](i, 0);
		ddGAp[i][1][0] = dataA[1][1](i, 0);
		ddGAp[i][1][1] = dataA[0][2](i, 0);
		dddGAp[i][0][0][0] = dataA[3][0](i, 0);
		dddGAp[i][0][0][1] = dataA[2][1](i, 0);
		dddGAp[i][0][1][0] = dataA[2][1](i, 0);
		dddGAp[i][0][1][1] = dataA[1][2](i, 0);
		dddGAp[i][1][0][0] = dataA[2][1](i, 0);
		dddGAp[i][1][0][1] = dataA[1][2](i, 0);
		dddGAp[i][1][1][0] = dataA[1][2](i, 0);
		dddGAp[i][1][1][1] = dataA[0][3](i, 0);

		GBp[i] = dataB[0][0](i, 0);
		dGBp[i][0] = dataB[1][0](i, 0);
		dGBp[i][1] = dataB[0][1](i, 0);
		ddGBp[i][0][0] = dataB[2][0](i, 0);
		ddGBp[i][0][1] = dataB[1][1](i, 0);
		ddGBp[i][1][0] = dataB[1][1](i, 0);
		ddGBp[i][1][1] = dataB[0][2](i, 0);
		dddGBp[i][0][0][0] = dataB[3][0](i, 0);
		dddGBp[i][0][0][1] = dataB[2][1](i, 0);
		dddGBp[i][0][1][0] = dataB[2][1](i, 0);
		dddGBp[i][0][1][1] = dataB[1][2](i, 0);
		dddGBp[i][1][0][0] = dataB[2][1](i, 0);
		dddGBp[i][1][0][1] = dataB[1][2](i, 0);
		dddGBp[i][1][1][0] = dataB[1][2](i, 0);
		dddGBp[i][1][1][1] = dataB[0][3](i, 0);
	}
}

//Evaluates the necessary derivatives to be used as input to AceGen routines
void RigidNURBS_1_RigidNURBS_1::EvaluateNURBSDerivatives_i(Matrix& mc)
{
	if (mc == last_ci && first_ci == false)
		return;
	else
	{
		last_ci = mc;
		first_ci = false;
	}
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	RigidNURBS_1* surfA;		//Ponteiro para a superfície A
	RigidNURBS_1* surfB;		//Ponteiro para a superfície B
	surfA = static_cast<RigidNURBS_1*>(db.surfaces[surf1_ID - 1]);
	surfB = static_cast<RigidNURBS_1*>(db.surfaces[surf2_ID - 1]);
	//Evaluates data and store in dataA and dataB
	int d = 1;
	surfA->surf->NURBSDerivatives(mc(0, 0), mc(1, 0), dataA, d);
	surfB->surf->NURBSDerivatives(mc(2, 0), mc(3, 0), dataB, d);
	//Copies data from dataA and dataB to AceGen input data pointers
	cAi[0] = mc(0, 0);
	cAi[1] = mc(1, 0);
	cBi[0] = mc(2, 0);
	cBi[1] = mc(3, 0);
	for (int i = 0; i < 3; i++)
	{
		GAi[i] = dataA[0][0](i, 0);
		dGAi[i][0] = dataA[1][0](i, 0);
		dGAi[i][1] = dataA[0][1](i, 0);
		
		GBi[i] = dataB[0][0](i, 0);
		dGBi[i][0] = dataB[1][0](i, 0);
		dGBi[i][1] = dataB[0][1](i, 0);
	}
}

//Evaluates the necessary DOFs variables to be used as input to AceGen routines
void RigidNURBS_1_RigidNURBS_1::EvaluateNURBSDOFsVariables()
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	RigidNURBS_1* surfA;		//Ponteiro para a superfície A
	RigidNURBS_1* surfB;		//Ponteiro para a superfície B
	surfA = static_cast<RigidNURBS_1*>(db.surfaces[surf1_ID - 1]);
	surfB = static_cast<RigidNURBS_1*>(db.surfaces[surf2_ID - 1]);

	for (int i = 0; i < 3; i++)
	{
		xAi[i] = surfA->xi(i, 0);
		xBi[i] = surfB->xi(i, 0);

		uA[i] = db.nodes[surfA->nodes[0] - 1]->displacements[i];
		uB[i] = db.nodes[surfB->nodes[0] - 1]->displacements[i];

		alphaA[i] = db.nodes[surfA->nodes[0] - 1]->displacements[i + 3];
		alphaB[i] = db.nodes[surfB->nodes[0] - 1]->displacements[i + 3];


		duiA[i] = db.nodes[surfA->nodes[0] - 1]->copy_vel[i];
		duiB[i] = db.nodes[surfB->nodes[0] - 1]->copy_vel[i];
		dduiA[i] = db.nodes[surfA->nodes[0] - 1]->copy_accel[i];
		dduiB[i] = db.nodes[surfB->nodes[0] - 1]->copy_accel[i];
		dalphaiA[i] = db.nodes[surfA->nodes[0] - 1]->copy_vel[i + 3];
		dalphaiB[i] = db.nodes[surfB->nodes[0] - 1]->copy_vel[i + 3];
		ddalphaiA[i] = db.nodes[surfA->nodes[0] - 1]->copy_accel[i + 3];
		ddalphaiB[i] = db.nodes[surfB->nodes[0] - 1]->copy_accel[i + 3];
		

		for (int j = 0; j < 3; j++)
		{
			QAi[i][j] = surfA->Qi(i, j);
			QBi[i][j] = surfB->Qi(i, j);
		}
	}
	invertnormalA = &surfA->invert_normal;
	invertnormalB = &surfB->invert_normal;
}

void RigidNURBS_1_RigidNURBS_1::ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft)
{
	//AceGen variables or pointers
	double v[15000];
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
	for (int i = 0; i < 12; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < 12; j++)
			Kc[i][j] = 0.0;
	}

	Matrix temp_cp(4);
	for (int i = 0; i < 4; i++)
		temp_cp(i, 0) = convective[i];
	Matrix temp_ci(4);
	for (int i = 0; i < 4; i++)
		temp_ci(i, 0) = copy_convective[i];

	EvaluateNURBSDerivatives_p(temp_cp);
	EvaluateNURBSDerivatives_i(temp_ci);
	EvaluateNURBSDOFsVariables();


#pragma region AceGen
	double v01; double v010; double v011; double v012; double v013; double v014;
	double v015; double v016; double v017; double v018; double v019; double v02;
	double v020; double v021; double v022; double v023; double v024; double v025;
	double v026; double v027; double v03; double v04; double v05; double v06; double v07;
	double v08; double v09;
	int i847, i915, i1143, i1429, i2302, i2362, i2773, i2774, i2775, i2776, i2777, i2778, i2819
		, i2820, i2821, i2822, i2823, i2824, i2833, i2834, i2835, i2845, i2846, i2847, i2848, i2849
		, i2850, i2992, i2993, i2994, i2995, i2996, i2997, b709, b710, b742, b776, b777, b778, b794
		, b811, b812, b824, b1092, b1148, b1149, b1150, b1191, b1228, b1687, b1732, b1800, b1801
		, b1818, b1936, b1961;
	v[2783] = -((*a4)*(*ct));
	v[4084] = 0e0;
	v[4085] = 0e0;
	v[4086] = 1e0;
	v[4087] = 0e0;
	v[4088] = 0e0;
	v[4089] = 0e0;
	v[4090] = 0e0;
	v[4091] = 0e0;
	v[4092] = -1e0;
	v[4093] = 0e0;
	v[4094] = 0e0;
	v[4095] = 0e0;
	v[4096] = 0e0;
	v[4097] = 1e0;
	v[4098] = 0e0;
	v[4099] = 0e0;
	v[4100] = 0e0;
	v[4101] = 0e0;
	v[4102] = 0e0;
	v[4103] = -1e0;
	v[4104] = 0e0;
	v[4105] = 0e0;
	v[4106] = 0e0;
	v[4107] = 0e0;
	v[4108] = 1e0;
	v[4109] = 0e0;
	v[4110] = 0e0;
	v[4111] = 0e0;
	v[4112] = 0e0;
	v[4113] = 0e0;
	v[4114] = -1e0;
	v[4115] = 0e0;
	v[4116] = 0e0;
	v[4117] = 0e0;
	v[4118] = 0e0;
	v[4119] = 0e0;
	v[868] = GAp[0] * QAi[0][0] + GAp[1] * QAi[0][1] + GAp[2] * QAi[0][2];
	v[866] = GAp[0] * QAi[1][0] + GAp[1] * QAi[1][1] + GAp[2] * QAi[1][2];
	v[864] = GAp[0] * QAi[2][0] + GAp[1] * QAi[2][1] + GAp[2] * QAi[2][2];
	v[401] = alphaA[0] / 2e0;
	v[399] = 2e0*alphaA[0];
	v[226] = Power(alphaA[0], 2);
	v[402] = 2e0*alphaA[1];
	v[3964] = 0e0;
	v[3965] = 0e0;
	v[3966] = 0e0;
	v[3967] = v[399];
	v[3968] = v[402];
	v[3969] = 0e0;
	v[3970] = 0e0;
	v[3971] = 0e0;
	v[3972] = 0e0;
	v[3973] = 0e0;
	v[3974] = 0e0;
	v[3975] = 0e0;
	v[4260] = 0e0;
	v[4261] = 0e0;
	v[4262] = 0e0;
	v[4263] = -v[399];
	v[4264] = -v[402];
	v[4265] = 0e0;
	v[4266] = 0e0;
	v[4267] = 0e0;
	v[4268] = 0e0;
	v[4269] = 0e0;
	v[4270] = 0e0;
	v[4271] = 0e0;
	v[400] = alphaA[1] / 2e0;
	v[3784] = 0e0;
	v[3785] = 0e0;
	v[3786] = 0e0;
	v[3787] = v[400];
	v[3788] = v[401];
	v[3789] = 0e0;
	v[3790] = 0e0;
	v[3791] = 0e0;
	v[3792] = 0e0;
	v[3793] = 0e0;
	v[3794] = 0e0;
	v[3795] = 0e0;
	v[224] = (alphaA[0] * alphaA[1]) / 2e0;
	v[219] = Power(alphaA[1], 2);
	v[451] = -v[219] - v[226];
	v[429] = alphaA[2] + v[224];
	v[419] = -alphaA[2] + v[224];
	v[404] = 2e0*alphaA[2];
	v[5120] = 0e0;
	v[5121] = 0e0;
	v[5122] = 0e0;
	v[5123] = -v[399];
	v[5124] = 0e0;
	v[5125] = -v[404];
	v[5126] = 0e0;
	v[5127] = 0e0;
	v[5128] = 0e0;
	v[5129] = 0e0;
	v[5130] = 0e0;
	v[5131] = 0e0;
	v[3916] = 0e0;
	v[3917] = 0e0;
	v[3918] = 0e0;
	v[3919] = v[399];
	v[3920] = 0e0;
	v[3921] = v[404];
	v[3922] = 0e0;
	v[3923] = 0e0;
	v[3924] = 0e0;
	v[3925] = 0e0;
	v[3926] = 0e0;
	v[3927] = 0e0;
	v[3868] = 0e0;
	v[3869] = 0e0;
	v[3870] = 0e0;
	v[3871] = 0e0;
	v[3872] = v[402];
	v[3873] = v[404];
	v[3874] = 0e0;
	v[3875] = 0e0;
	v[3876] = 0e0;
	v[3877] = 0e0;
	v[3878] = 0e0;
	v[3879] = 0e0;
	v[4236] = 0e0;
	v[4237] = 0e0;
	v[4238] = 0e0;
	v[4239] = 0e0;
	v[4240] = -v[402];
	v[4241] = -v[404];
	v[4242] = 0e0;
	v[4243] = 0e0;
	v[4244] = 0e0;
	v[4245] = 0e0;
	v[4246] = 0e0;
	v[4247] = 0e0;
	v[3772] = 0e0;
	v[3773] = 0e0;
	v[3774] = 0e0;
	v[3775] = v[399];
	v[3776] = v[402];
	v[3777] = v[404];
	v[3778] = 0e0;
	v[3779] = 0e0;
	v[3780] = 0e0;
	v[3781] = 0e0;
	v[3782] = 0e0;
	v[3783] = 0e0;
	v[403] = alphaA[2] / 2e0;
	v[3796] = 0e0;
	v[3797] = 0e0;
	v[3798] = 0e0;
	v[3799] = v[403];
	v[3800] = 0e0;
	v[3801] = v[401];
	v[3802] = 0e0;
	v[3803] = 0e0;
	v[3804] = 0e0;
	v[3805] = 0e0;
	v[3806] = 0e0;
	v[3807] = 0e0;
	v[3808] = 0e0;
	v[3809] = 0e0;
	v[3810] = 0e0;
	v[3811] = 0e0;
	v[3812] = v[403];
	v[3813] = v[400];
	v[3814] = 0e0;
	v[3815] = 0e0;
	v[3816] = 0e0;
	v[3817] = 0e0;
	v[3818] = 0e0;
	v[3819] = 0e0;
	v[231] = (alphaA[1] * alphaA[2]) / 2e0;
	v[446] = alphaA[0] + v[231];
	v[438] = -alphaA[0] + v[231];
	v[229] = (alphaA[0] * alphaA[2]) / 2e0;
	v[442] = -alphaA[1] + v[229];
	v[425] = alphaA[1] + v[229];
	v[220] = Power(alphaA[2], 2);
	v[1068] = 4e0 + v[219] + v[220] + v[226];
	v[3173] = 1e0 / Power(v[1068], 4);
	v[1411] = 1e0 / Power(v[1068], 3);
	v[1484] = -8e0*v[1411] * v[399];
	v[1482] = -8e0*v[1411] * v[404];
	v[1480] = 8e0*v[1411] * v[402];
	v[962] = 1e0 / Power(v[1068], 2);
	v[2959] = 4e0*v[962];
	v[433] = -v[220] - v[226];
	v[414] = -v[219] - v[220];
	v[2410] = -2e0*v[414] * v[962];
	v[413] = 4e0*v[404] * v[962];
	v[455] = -(v[413] * v[451]) / 2e0;
	v[412] = -4e0*v[402] * v[962];
	v[435] = (v[412] * v[433]) / 2e0;
	v[411] = 4e0*v[399] * v[962];
	v[415] = -(v[411] * v[414]) / 2e0;
	v[310] = (*a4)*alphaA[0] + (*a5)*dalphaiA[0] + (*a6)*ddalphaiA[0];
	v[316] = (*a4)*alphaA[1] + (*a5)*dalphaiA[1] + (*a6)*ddalphaiA[1];
	v[318] = (*a4)*alphaA[2] + (*a5)*dalphaiA[2] + (*a6)*ddalphaiA[2];
	v[856] = -(GBp[0] * QBi[0][0]) - GBp[1] * QBi[0][1] - GBp[2] * QBi[0][2];
	v[854] = -(GBp[0] * QBi[1][0]) - GBp[1] * QBi[1][1] - GBp[2] * QBi[1][2];
	v[852] = -(GBp[0] * QBi[2][0]) - GBp[1] * QBi[2][1] - GBp[2] * QBi[2][2];
	v[407] = alphaB[0] / 2e0;
	v[405] = 2e0*alphaB[0];
	v[266] = Power(alphaB[0], 2);
	v[408] = 2e0*alphaB[1];
	v[4072] = 0e0;
	v[4073] = 0e0;
	v[4074] = 0e0;
	v[4075] = 0e0;
	v[4076] = 0e0;
	v[4077] = 0e0;
	v[4078] = 0e0;
	v[4079] = 0e0;
	v[4080] = 0e0;
	v[4081] = v[405];
	v[4082] = v[408];
	v[4083] = 0e0;
	v[4332] = 0e0;
	v[4333] = 0e0;
	v[4334] = 0e0;
	v[4335] = 0e0;
	v[4336] = 0e0;
	v[4337] = 0e0;
	v[4338] = 0e0;
	v[4339] = 0e0;
	v[4340] = 0e0;
	v[4341] = -v[405];
	v[4342] = -v[408];
	v[4343] = 0e0;
	v[406] = alphaB[1] / 2e0;
	v[3832] = 0e0;
	v[3833] = 0e0;
	v[3834] = 0e0;
	v[3835] = 0e0;
	v[3836] = 0e0;
	v[3837] = 0e0;
	v[3838] = 0e0;
	v[3839] = 0e0;
	v[3840] = 0e0;
	v[3841] = v[406];
	v[3842] = v[407];
	v[3843] = 0e0;
	v[264] = (alphaB[0] * alphaB[1]) / 2e0;
	v[259] = Power(alphaB[1], 2);
	v[550] = -v[259] - v[266];
	v[528] = alphaB[2] + v[264];
	v[518] = -alphaB[2] + v[264];
	v[410] = 2e0*alphaB[2];
	v[5132] = 0e0;
	v[5133] = 0e0;
	v[5134] = 0e0;
	v[5135] = 0e0;
	v[5136] = 0e0;
	v[5137] = 0e0;
	v[5138] = 0e0;
	v[5139] = 0e0;
	v[5140] = 0e0;
	v[5141] = -v[405];
	v[5142] = 0e0;
	v[5143] = -v[410];
	v[4024] = 0e0;
	v[4025] = 0e0;
	v[4026] = 0e0;
	v[4027] = 0e0;
	v[4028] = 0e0;
	v[4029] = 0e0;
	v[4030] = 0e0;
	v[4031] = 0e0;
	v[4032] = 0e0;
	v[4033] = v[405];
	v[4034] = 0e0;
	v[4035] = v[410];
	v[3976] = 0e0;
	v[3977] = 0e0;
	v[3978] = 0e0;
	v[3979] = 0e0;
	v[3980] = 0e0;
	v[3981] = 0e0;
	v[3982] = 0e0;
	v[3983] = 0e0;
	v[3984] = 0e0;
	v[3985] = 0e0;
	v[3986] = v[408];
	v[3987] = v[410];
	v[4308] = 0e0;
	v[4309] = 0e0;
	v[4310] = 0e0;
	v[4311] = 0e0;
	v[4312] = 0e0;
	v[4313] = 0e0;
	v[4314] = 0e0;
	v[4315] = 0e0;
	v[4316] = 0e0;
	v[4317] = 0e0;
	v[4318] = -v[408];
	v[4319] = -v[410];
	v[3820] = 0e0;
	v[3821] = 0e0;
	v[3822] = 0e0;
	v[3823] = 0e0;
	v[3824] = 0e0;
	v[3825] = 0e0;
	v[3826] = 0e0;
	v[3827] = 0e0;
	v[3828] = 0e0;
	v[3829] = v[405];
	v[3830] = v[408];
	v[3831] = v[410];
	v[409] = alphaB[2] / 2e0;
	v[3844] = 0e0;
	v[3845] = 0e0;
	v[3846] = 0e0;
	v[3847] = 0e0;
	v[3848] = 0e0;
	v[3849] = 0e0;
	v[3850] = 0e0;
	v[3851] = 0e0;
	v[3852] = 0e0;
	v[3853] = v[409];
	v[3854] = 0e0;
	v[3855] = v[407];
	v[3856] = 0e0;
	v[3857] = 0e0;
	v[3858] = 0e0;
	v[3859] = 0e0;
	v[3860] = 0e0;
	v[3861] = 0e0;
	v[3862] = 0e0;
	v[3863] = 0e0;
	v[3864] = 0e0;
	v[3865] = 0e0;
	v[3866] = v[409];
	v[3867] = v[406];
	v[271] = (alphaB[1] * alphaB[2]) / 2e0;
	v[545] = alphaB[0] + v[271];
	v[537] = -alphaB[0] + v[271];
	v[269] = (alphaB[0] * alphaB[2]) / 2e0;
	v[541] = -alphaB[1] + v[269];
	v[524] = alphaB[1] + v[269];
	v[260] = Power(alphaB[2], 2);
	v[1078] = 4e0 + v[259] + v[260] + v[266];
	v[3167] = 1e0 / Power(v[1078], 4);
	v[1413] = 1e0 / Power(v[1078], 3);
	v[1535] = -8e0*v[1413] * v[405];
	v[1533] = -8e0*v[1413] * v[410];
	v[1531] = 8e0*v[1413] * v[408];
	v[977] = 1e0 / Power(v[1078], 2);
	v[2963] = 4e0*v[977];
	v[532] = -v[260] - v[266];
	v[513] = -v[259] - v[260];
	v[2421] = -2e0*v[513] * v[977];
	v[512] = 4e0*v[410] * v[977];
	v[554] = -(v[512] * v[550]) / 2e0;
	v[511] = -4e0*v[408] * v[977];
	v[534] = (v[511] * v[532]) / 2e0;
	v[510] = 4e0*v[405] * v[977];
	v[514] = -(v[510] * v[513]) / 2e0;
	v[336] = (*a4)*alphaB[0] + (*a5)*dalphaiB[0] + (*a6)*ddalphaiB[0];
	v[342] = (*a4)*alphaB[1] + (*a5)*dalphaiB[1] + (*a6)*ddalphaiB[1];
	v[344] = (*a4)*alphaB[2] + (*a5)*dalphaiB[2] + (*a6)*ddalphaiB[2];
	v[218] = 4e0 / v[1068];
	v[2957] = -v[218] / 2e0;
	v[453] = -(v[218] * v[402]) / 2e0;
	v[454] = (v[412] * v[451]) / 2e0 + v[453];
	v[450] = -(v[218] * v[399]) / 2e0;
	v[452] = v[450] - (v[411] * v[451]) / 2e0;
	v[447] = v[218] - v[411] * v[446];
	v[444] = -v[218] + v[412] * v[442];
	v[439] = -v[218] - v[411] * v[438];
	v[436] = -(v[218] * v[404]) / 2e0;
	v[437] = -(v[413] * v[433]) / 2e0 + v[436];
	v[434] = -(v[411] * v[433]) / 2e0 + v[450];
	v[432] = v[218] - v[413] * v[429];
	v[427] = v[218] + v[412] * v[425];
	v[424] = -(v[218] * v[403]);
	v[448] = -v[424] + v[412] * v[446];
	v[493] = QAi[0][2] * v[444] + QAi[1][2] * v[448] + QAi[2][2] * v[454];
	v[490] = QAi[0][1] * v[444] + QAi[1][1] * v[448] + QAi[2][1] * v[454];
	v[487] = QAi[0][0] * v[444] + QAi[1][0] * v[448] + QAi[2][0] * v[454];
	v[508] = GAp[0] * v[487] + GAp[1] * v[490] + GAp[2] * v[493];
	v[443] = -v[424] - v[411] * v[442];
	v[492] = QAi[0][2] * v[443] + QAi[1][2] * v[447] + QAi[2][2] * v[452];
	v[489] = QAi[0][1] * v[443] + QAi[1][1] * v[447] + QAi[2][1] * v[452];
	v[486] = QAi[0][0] * v[443] + QAi[1][0] * v[447] + QAi[2][0] * v[452];
	v[507] = GAp[0] * v[486] + GAp[1] * v[489] + GAp[2] * v[492];
	v[440] = -v[424] + v[412] * v[438];
	v[426] = -v[424] - v[411] * v[425];
	v[423] = -v[218] - v[413] * v[419];
	v[421] = -(v[218] * v[401]);
	v[445] = -v[421] - v[413] * v[442];
	v[431] = -v[421] + v[412] * v[429];
	v[478] = QAi[0][2] * v[431] + QAi[1][2] * v[435] + QAi[2][2] * v[440];
	v[475] = QAi[0][1] * v[431] + QAi[1][1] * v[435] + QAi[2][1] * v[440];
	v[472] = QAi[0][0] * v[431] + QAi[1][0] * v[435] + QAi[2][0] * v[440];
	v[505] = GAp[0] * v[472] + GAp[1] * v[475] + GAp[2] * v[478];
	v[428] = -v[421] - v[413] * v[425];
	v[422] = v[412] * v[419] - v[421];
	v[418] = v[218] * v[400];
	v[449] = v[418] - v[413] * v[446];
	v[494] = QAi[0][2] * v[445] + QAi[1][2] * v[449] + QAi[2][2] * v[455];
	v[491] = QAi[0][1] * v[445] + QAi[1][1] * v[449] + QAi[2][1] * v[455];
	v[488] = QAi[0][0] * v[445] + QAi[1][0] * v[449] + QAi[2][0] * v[455];
	v[509] = GAp[0] * v[488] + GAp[1] * v[491] + GAp[2] * v[494];
	v[441] = v[418] - v[413] * v[438];
	v[479] = QAi[0][2] * v[432] + QAi[1][2] * v[437] + QAi[2][2] * v[441];
	v[476] = QAi[0][1] * v[432] + QAi[1][1] * v[437] + QAi[2][1] * v[441];
	v[473] = QAi[0][0] * v[432] + QAi[1][0] * v[437] + QAi[2][0] * v[441];
	v[506] = GAp[0] * v[473] + GAp[1] * v[476] + GAp[2] * v[479];
	v[430] = v[418] - v[411] * v[429];
	v[477] = QAi[0][2] * v[430] + QAi[1][2] * v[434] + QAi[2][2] * v[439];
	v[474] = QAi[0][1] * v[430] + QAi[1][1] * v[434] + QAi[2][1] * v[439];
	v[471] = QAi[0][0] * v[430] + QAi[1][0] * v[434] + QAi[2][0] * v[439];
	v[504] = GAp[0] * v[471] + GAp[1] * v[474] + GAp[2] * v[477];
	v[420] = v[418] - v[411] * v[419];
	v[462] = QAi[0][2] * v[415] + QAi[1][2] * v[420] + QAi[2][2] * v[426];
	v[459] = QAi[0][1] * v[415] + QAi[1][1] * v[420] + QAi[2][1] * v[426];
	v[456] = QAi[0][0] * v[415] + QAi[1][0] * v[420] + QAi[2][0] * v[426];
	v[501] = GAp[0] * v[456] + GAp[1] * v[459] + GAp[2] * v[462];
	v[417] = -(v[413] * v[414]) / 2e0 + v[436];
	v[464] = QAi[0][2] * v[417] + QAi[1][2] * v[423] + QAi[2][2] * v[428];
	v[461] = QAi[0][1] * v[417] + QAi[1][1] * v[423] + QAi[2][1] * v[428];
	v[458] = QAi[0][0] * v[417] + QAi[1][0] * v[423] + QAi[2][0] * v[428];
	v[503] = GAp[0] * v[458] + GAp[1] * v[461] + GAp[2] * v[464];
	v[416] = (v[412] * v[414]) / 2e0 + v[453];
	v[463] = QAi[0][2] * v[416] + QAi[1][2] * v[422] + QAi[2][2] * v[427];
	v[460] = QAi[0][1] * v[416] + QAi[1][1] * v[422] + QAi[2][1] * v[427];
	v[457] = QAi[0][0] * v[416] + QAi[1][0] * v[422] + QAi[2][0] * v[427];
	v[502] = GAp[0] * v[457] + GAp[1] * v[460] + GAp[2] * v[463];
	v[307] = (v[218] * v[218]);
	v[3177] = 1e0 / Power(v[307], 3);
	v[2705] = v[318] / v[307];
	v[2704] = v[316] / v[307];
	v[2703] = v[310] / v[307];
	v[221] = 1e0 + (v[218] * v[414]) / 2e0;
	v[1314] = v[221] * v[2703];
	v[222] = v[218] * v[419];
	v[1316] = v[222] * v[2704];
	v[1764] = v[1314] + v[1316];
	v[223] = v[218] * v[425];
	v[3002] = v[223] / v[307];
	v[1317] = v[223] * v[2705];
	v[1769] = v[1314] + v[1317];
	v[1758] = v[1316] + v[1769];
	v[225] = v[218] * v[429];
	v[1321] = v[225] * v[2703];
	v[227] = 1e0 + (v[218] * v[433]) / 2e0;
	v[1309] = v[227] * v[2704];
	v[1760] = v[1309] + v[1321];
	v[228] = v[218] * v[438];
	v[3003] = v[228] / v[307];
	v[1310] = v[228] * v[2705];
	v[1770] = v[1309] + v[1310];
	v[1763] = v[1321] + v[1770];
	v[230] = v[218] * v[442];
	v[1325] = v[230] * v[2703];
	v[232] = v[218] * v[446];
	v[3004] = v[232] / v[307];
	v[1305] = v[232] * v[2704];
	v[233] = 1e0 + (v[218] * v[451]) / 2e0;
	v[1306] = v[233] * v[2705];
	v[1768] = v[1305] + v[1306] + v[1325];
	v[1765] = -v[1325] + v[1768];
	v[1759] = -v[1305] + v[1768];
	v[325] = -(v[418] * v[424]);
	v[2713] = v[325] - v[411];
	v[323] = v[421] * v[424];
	v[2711] = v[323] - v[412];
	v[314] = -(v[418] * v[421]);
	v[2709] = v[314] - v[413];
	v[237] = QAi[0][0] * v[221] + QAi[1][0] * v[222] + QAi[2][0] * v[223];
	v[238] = QAi[0][1] * v[221] + QAi[1][1] * v[222] + QAi[2][1] * v[223];
	v[239] = QAi[0][2] * v[221] + QAi[1][2] * v[222] + QAi[2][2] * v[223];
	v[387] = dGAp[0][1] * v[237] + dGAp[1][1] * v[238] + dGAp[2][1] * v[239];
	v[383] = dGAp[0][0] * v[237] + dGAp[1][0] * v[238] + dGAp[2][0] * v[239];
	v[240] = QAi[0][0] * v[225] + QAi[1][0] * v[227] + QAi[2][0] * v[228];
	v[241] = QAi[0][1] * v[225] + QAi[1][1] * v[227] + QAi[2][1] * v[228];
	v[242] = QAi[0][2] * v[225] + QAi[1][2] * v[227] + QAi[2][2] * v[228];
	v[388] = dGAp[0][1] * v[240] + dGAp[1][1] * v[241] + dGAp[2][1] * v[242];
	v[384] = dGAp[0][0] * v[240] + dGAp[1][0] * v[241] + dGAp[2][0] * v[242];
	v[243] = QAi[0][0] * v[230] + QAi[1][0] * v[232] + QAi[2][0] * v[233];
	v[244] = QAi[0][1] * v[230] + QAi[1][1] * v[232] + QAi[2][1] * v[233];
	v[245] = QAi[0][2] * v[230] + QAi[1][2] * v[232] + QAi[2][2] * v[233];
	v[389] = dGAp[0][1] * v[243] + dGAp[1][1] * v[244] + dGAp[2][1] * v[245];
	v[385] = dGAp[0][0] * v[243] + dGAp[1][0] * v[244] + dGAp[2][0] * v[245];
	v[258] = 4e0 / v[1078];
	v[2961] = -v[258] / 2e0;
	v[552] = -(v[258] * v[408]) / 2e0;
	v[553] = (v[511] * v[550]) / 2e0 + v[552];
	v[549] = -(v[258] * v[405]) / 2e0;
	v[551] = v[549] - (v[510] * v[550]) / 2e0;
	v[546] = v[258] - v[510] * v[545];
	v[543] = -v[258] + v[511] * v[541];
	v[538] = -v[258] - v[510] * v[537];
	v[535] = -(v[258] * v[410]) / 2e0;
	v[536] = -(v[512] * v[532]) / 2e0 + v[535];
	v[533] = -(v[510] * v[532]) / 2e0 + v[549];
	v[531] = v[258] - v[512] * v[528];
	v[526] = v[258] + v[511] * v[524];
	v[523] = -(v[258] * v[409]);
	v[547] = -v[523] + v[511] * v[545];
	v[592] = QBi[0][2] * v[543] + QBi[1][2] * v[547] + QBi[2][2] * v[553];
	v[589] = QBi[0][1] * v[543] + QBi[1][1] * v[547] + QBi[2][1] * v[553];
	v[586] = QBi[0][0] * v[543] + QBi[1][0] * v[547] + QBi[2][0] * v[553];
	v[607] = GBp[0] * v[586] + GBp[1] * v[589] + GBp[2] * v[592];
	v[542] = -v[523] - v[510] * v[541];
	v[591] = QBi[0][2] * v[542] + QBi[1][2] * v[546] + QBi[2][2] * v[551];
	v[588] = QBi[0][1] * v[542] + QBi[1][1] * v[546] + QBi[2][1] * v[551];
	v[585] = QBi[0][0] * v[542] + QBi[1][0] * v[546] + QBi[2][0] * v[551];
	v[606] = GBp[0] * v[585] + GBp[1] * v[588] + GBp[2] * v[591];
	v[539] = -v[523] + v[511] * v[537];
	v[525] = -v[523] - v[510] * v[524];
	v[522] = -v[258] - v[512] * v[518];
	v[520] = -(v[258] * v[407]);
	v[544] = -v[520] - v[512] * v[541];
	v[530] = -v[520] + v[511] * v[528];
	v[577] = QBi[0][2] * v[530] + QBi[1][2] * v[534] + QBi[2][2] * v[539];
	v[574] = QBi[0][1] * v[530] + QBi[1][1] * v[534] + QBi[2][1] * v[539];
	v[571] = QBi[0][0] * v[530] + QBi[1][0] * v[534] + QBi[2][0] * v[539];
	v[604] = GBp[0] * v[571] + GBp[1] * v[574] + GBp[2] * v[577];
	v[527] = -v[520] - v[512] * v[524];
	v[521] = v[511] * v[518] - v[520];
	v[517] = v[258] * v[406];
	v[548] = v[517] - v[512] * v[545];
	v[593] = QBi[0][2] * v[544] + QBi[1][2] * v[548] + QBi[2][2] * v[554];
	v[590] = QBi[0][1] * v[544] + QBi[1][1] * v[548] + QBi[2][1] * v[554];
	v[587] = QBi[0][0] * v[544] + QBi[1][0] * v[548] + QBi[2][0] * v[554];
	v[608] = GBp[0] * v[587] + GBp[1] * v[590] + GBp[2] * v[593];
	v[540] = v[517] - v[512] * v[537];
	v[578] = QBi[0][2] * v[531] + QBi[1][2] * v[536] + QBi[2][2] * v[540];
	v[575] = QBi[0][1] * v[531] + QBi[1][1] * v[536] + QBi[2][1] * v[540];
	v[572] = QBi[0][0] * v[531] + QBi[1][0] * v[536] + QBi[2][0] * v[540];
	v[605] = GBp[0] * v[572] + GBp[1] * v[575] + GBp[2] * v[578];
	v[529] = v[517] - v[510] * v[528];
	v[576] = QBi[0][2] * v[529] + QBi[1][2] * v[533] + QBi[2][2] * v[538];
	v[573] = QBi[0][1] * v[529] + QBi[1][1] * v[533] + QBi[2][1] * v[538];
	v[570] = QBi[0][0] * v[529] + QBi[1][0] * v[533] + QBi[2][0] * v[538];
	v[603] = GBp[0] * v[570] + GBp[1] * v[573] + GBp[2] * v[576];
	v[519] = v[517] - v[510] * v[518];
	v[561] = QBi[0][2] * v[514] + QBi[1][2] * v[519] + QBi[2][2] * v[525];
	v[558] = QBi[0][1] * v[514] + QBi[1][1] * v[519] + QBi[2][1] * v[525];
	v[555] = QBi[0][0] * v[514] + QBi[1][0] * v[519] + QBi[2][0] * v[525];
	v[600] = GBp[0] * v[555] + GBp[1] * v[558] + GBp[2] * v[561];
	v[618] = -(v[387] * v[600]) - v[388] * v[603] - v[389] * v[606];
	v[612] = -(v[383] * v[600]) - v[384] * v[603] - v[385] * v[606];
	v[516] = -(v[512] * v[513]) / 2e0 + v[535];
	v[563] = QBi[0][2] * v[516] + QBi[1][2] * v[522] + QBi[2][2] * v[527];
	v[560] = QBi[0][1] * v[516] + QBi[1][1] * v[522] + QBi[2][1] * v[527];
	v[557] = QBi[0][0] * v[516] + QBi[1][0] * v[522] + QBi[2][0] * v[527];
	v[602] = GBp[0] * v[557] + GBp[1] * v[560] + GBp[2] * v[563];
	v[620] = -(v[387] * v[602]) - v[388] * v[605] - v[389] * v[608];
	v[614] = -(v[383] * v[602]) - v[384] * v[605] - v[385] * v[608];
	v[515] = (v[511] * v[513]) / 2e0 + v[552];
	v[562] = QBi[0][2] * v[515] + QBi[1][2] * v[521] + QBi[2][2] * v[526];
	v[559] = QBi[0][1] * v[515] + QBi[1][1] * v[521] + QBi[2][1] * v[526];
	v[556] = QBi[0][0] * v[515] + QBi[1][0] * v[521] + QBi[2][0] * v[526];
	v[601] = GBp[0] * v[556] + GBp[1] * v[559] + GBp[2] * v[562];
	v[619] = -(v[387] * v[601]) - v[388] * v[604] - v[389] * v[607];
	v[613] = -(v[383] * v[601]) - v[384] * v[604] - v[385] * v[607];
	v[333] = (v[258] * v[258]);
	v[3171] = 1e0 / Power(v[333], 3);
	v[2708] = v[344] / v[333];
	v[2707] = v[342] / v[333];
	v[2706] = v[336] / v[333];
	v[261] = 1e0 + (v[258] * v[513]) / 2e0;
	v[1287] = v[261] * v[2706];
	v[262] = v[258] * v[518];
	v[1289] = v[262] * v[2707];
	v[1779] = v[1287] + v[1289];
	v[263] = v[258] * v[524];
	v[2999] = v[263] / v[333];
	v[1290] = v[263] * v[2708];
	v[1784] = v[1287] + v[1290];
	v[1773] = v[1289] + v[1784];
	v[265] = v[258] * v[528];
	v[1294] = v[265] * v[2706];
	v[267] = 1e0 + (v[258] * v[532]) / 2e0;
	v[1282] = v[267] * v[2707];
	v[1775] = v[1282] + v[1294];
	v[268] = v[258] * v[537];
	v[3000] = v[268] / v[333];
	v[1283] = v[268] * v[2708];
	v[1785] = v[1282] + v[1283];
	v[1778] = v[1294] + v[1785];
	v[270] = v[258] * v[541];
	v[1298] = v[270] * v[2706];
	v[272] = v[258] * v[545];
	v[3001] = v[272] / v[333];
	v[1278] = v[2707] * v[272];
	v[273] = 1e0 + (v[258] * v[550]) / 2e0;
	v[1279] = v[2708] * v[273];
	v[1783] = v[1278] + v[1279] + v[1298];
	v[1780] = -v[1298] + v[1783];
	v[1774] = -v[1278] + v[1783];
	v[351] = -(v[517] * v[523]);
	v[2718] = v[351] - v[510];
	v[349] = v[520] * v[523];
	v[2716] = v[349] - v[511];
	v[340] = -(v[517] * v[520]);
	v[2714] = v[340] - v[512];
	v[277] = QBi[0][0] * v[261] + QBi[1][0] * v[262] + QBi[2][0] * v[263];
	v[278] = QBi[0][1] * v[261] + QBi[1][1] * v[262] + QBi[2][1] * v[263];
	v[279] = QBi[0][2] * v[261] + QBi[1][2] * v[262] + QBi[2][2] * v[263];
	v[395] = dGBp[0][1] * v[277] + dGBp[1][1] * v[278] + dGBp[2][1] * v[279];
	v[391] = dGBp[0][0] * v[277] + dGBp[1][0] * v[278] + dGBp[2][0] * v[279];
	v[700] = -(invH[3][0] * v[383]) - invH[3][1] * v[387] + invH[3][2] * v[391] + invH[3][3] * v[395];
	v[685] = -(invH[2][0] * v[383]) - invH[2][1] * v[387] + invH[2][2] * v[391] + invH[2][3] * v[395];
	v[670] = -(invH[1][0] * v[383]) - invH[1][1] * v[387] + invH[1][2] * v[391] + invH[1][3] * v[395];
	v[655] = -(invH[0][0] * v[383]) - invH[0][1] * v[387] + invH[0][2] * v[391] + invH[0][3] * v[395];
	v[280] = QBi[0][0] * v[265] + QBi[1][0] * v[267] + QBi[2][0] * v[268];
	v[281] = QBi[0][1] * v[265] + QBi[1][1] * v[267] + QBi[2][1] * v[268];
	v[282] = QBi[0][2] * v[265] + QBi[1][2] * v[267] + QBi[2][2] * v[268];
	v[396] = dGBp[0][1] * v[280] + dGBp[1][1] * v[281] + dGBp[2][1] * v[282];
	v[392] = dGBp[0][0] * v[280] + dGBp[1][0] * v[281] + dGBp[2][0] * v[282];
	v[702] = -(invH[3][0] * v[384]) - invH[3][1] * v[388] + invH[3][2] * v[392] + invH[3][3] * v[396];
	v[687] = -(invH[2][0] * v[384]) - invH[2][1] * v[388] + invH[2][2] * v[392] + invH[2][3] * v[396];
	v[672] = -(invH[1][0] * v[384]) - invH[1][1] * v[388] + invH[1][2] * v[392] + invH[1][3] * v[396];
	v[657] = -(invH[0][0] * v[384]) - invH[0][1] * v[388] + invH[0][2] * v[392] + invH[0][3] * v[396];
	v[283] = QBi[0][0] * v[270] + QBi[1][0] * v[272] + QBi[2][0] * v[273];
	v[284] = QBi[0][1] * v[270] + QBi[1][1] * v[272] + QBi[2][1] * v[273];
	v[285] = QBi[0][2] * v[270] + QBi[1][2] * v[272] + QBi[2][2] * v[273];
	v[397] = dGBp[0][1] * v[283] + dGBp[1][1] * v[284] + dGBp[2][1] * v[285];
	v[629] = -(v[395] * v[503]) - v[396] * v[506] - v[397] * v[509];
	v[628] = -(v[395] * v[502]) - v[396] * v[505] - v[397] * v[508];
	v[627] = -(v[395] * v[501]) - v[396] * v[504] - v[397] * v[507];
	v[393] = dGBp[0][0] * v[283] + dGBp[1][0] * v[284] + dGBp[2][0] * v[285];
	v[704] = -(invH[3][0] * v[385]) - invH[3][1] * v[389] + invH[3][2] * v[393] + invH[3][3] * v[397];
	v[689] = -(invH[2][0] * v[385]) - invH[2][1] * v[389] + invH[2][2] * v[393] + invH[2][3] * v[397];
	v[674] = -(invH[1][0] * v[385]) - invH[1][1] * v[389] + invH[1][2] * v[393] + invH[1][3] * v[397];
	v[659] = -(invH[0][0] * v[385]) - invH[0][1] * v[389] + invH[0][2] * v[393] + invH[0][3] * v[397];
	v[623] = -(v[391] * v[503]) - v[392] * v[506] - v[393] * v[509];
	v[622] = -(v[391] * v[502]) - v[392] * v[505] - v[393] * v[508];
	v[621] = -(v[391] * v[501]) - v[392] * v[504] - v[393] * v[507];
	v[2765] = uA[0] - uB[0] + xAi[0] - xBi[0];
	v[2764] = uA[1] - uB[1] + xAi[1] - xBi[1];
	v[2763] = uA[2] - uB[2] + xAi[2] - xBi[2];
	v[298] = (*a6)*dduiA[0] + (*a5)*duiA[0] + (*a4)*uA[0];
	v[299] = (*a6)*dduiA[1] + (*a5)*duiA[1] + (*a4)*uA[1];
	v[300] = (*a6)*dduiA[2] + (*a5)*duiA[2] + (*a4)*uA[2];
	v[301] = (*a6)*dduiB[0] + (*a5)*duiB[0] + (*a4)*uB[0];
	v[2759] = -v[298] + v[301];
	v[302] = (*a6)*dduiB[1] + (*a5)*duiB[1] + (*a4)*uB[1];
	v[2760] = -v[299] + v[302];
	v[303] = (*a6)*dduiB[2] + (*a5)*duiB[2] + (*a4)*uB[2];
	v[2761] = -v[300] + v[303];
	v[305] = v[323] + v[412];
	v[2937] = v[233] * v[305];
	v[2712] = v[305] / v[307];
	v[306] = v[314] + v[413];
	v[2933] = v[227] * v[306];
	v[2710] = v[306] / v[307];
	v[308] = v[307] + (v[421] * v[421]);
	v[2293] = v[308] / v[307];
	v[1420] = v[2711] / v[307];
	v[1418] = v[2709] / v[307];
	v[4524] = 0e0;
	v[4525] = 0e0;
	v[4526] = 0e0;
	v[4527] = 0e0;
	v[4528] = v[1418];
	v[4529] = v[1420];
	v[4530] = 0e0;
	v[4531] = 0e0;
	v[4532] = 0e0;
	v[4533] = 0e0;
	v[4534] = 0e0;
	v[4535] = 0e0;
	v[309] = v[2293] * v[310] + v[1418] * v[316] + v[1420] * v[318];
	v[311] = v[307] + (v[418] * v[418]);
	v[2936] = v[228] * v[311];
	v[2935] = v[225] * v[311];
	v[2931] = v[221] * v[2709];
	v[317] = v[325] + v[411];
	v[2938] = v[233] * v[317];
	v[2294] = v[311] / v[307];
	v[1419] = v[2713] / v[307];
	v[4536] = 0e0;
	v[4537] = 0e0;
	v[4538] = 0e0;
	v[4539] = v[2710];
	v[4540] = 0e0;
	v[4541] = v[1419];
	v[4542] = 0e0;
	v[4543] = 0e0;
	v[4544] = 0e0;
	v[4545] = 0e0;
	v[4546] = 0e0;
	v[4547] = 0e0;
	v[319] = v[2710] * v[310] + v[2294] * v[316] + v[1419] * v[318];
	v[320] = v[307] + (v[424] * v[424]);
	v[2930] = v[232] * v[320];
	v[2929] = v[230] * v[320];
	v[2932] = v[221] * v[2711];
	v[2934] = v[227] * v[2713];
	v[2295] = v[320] / v[307];
	v[1417] = v[317] / v[307];
	v[4548] = 0e0;
	v[4549] = 0e0;
	v[4550] = 0e0;
	v[4551] = v[2712];
	v[4552] = v[1417];
	v[4553] = 0e0;
	v[4554] = 0e0;
	v[4555] = 0e0;
	v[4556] = 0e0;
	v[4557] = 0e0;
	v[4558] = 0e0;
	v[4559] = 0e0;
	v[329] = v[2712] * v[310] + v[1417] * v[316] + v[2295] * v[318];
	v[331] = v[349] + v[511];
	v[2924] = v[273] * v[331];
	v[2717] = v[331] / v[333];
	v[332] = v[340] + v[512];
	v[2920] = v[267] * v[332];
	v[2715] = v[332] / v[333];
	v[334] = v[333] + (v[520] * v[520]);
	v[2296] = v[334] / v[333];
	v[1426] = v[2716] / v[333];
	v[1424] = v[2714] / v[333];
	v[4560] = 0e0;
	v[4561] = 0e0;
	v[4562] = 0e0;
	v[4563] = 0e0;
	v[4564] = 0e0;
	v[4565] = 0e0;
	v[4566] = 0e0;
	v[4567] = 0e0;
	v[4568] = 0e0;
	v[4569] = 0e0;
	v[4570] = v[1424];
	v[4571] = v[1426];
	v[335] = v[2296] * v[336] + v[1424] * v[342] + v[1426] * v[344];
	v[337] = v[333] + (v[517] * v[517]);
	v[2923] = v[268] * v[337];
	v[2922] = v[265] * v[337];
	v[2918] = v[261] * v[2714];
	v[343] = v[351] + v[510];
	v[2925] = v[273] * v[343];
	v[2297] = v[337] / v[333];
	v[1425] = v[2718] / v[333];
	v[4572] = 0e0;
	v[4573] = 0e0;
	v[4574] = 0e0;
	v[4575] = 0e0;
	v[4576] = 0e0;
	v[4577] = 0e0;
	v[4578] = 0e0;
	v[4579] = 0e0;
	v[4580] = 0e0;
	v[4581] = v[2715];
	v[4582] = 0e0;
	v[4583] = v[1425];
	v[345] = v[2715] * v[336] + v[2297] * v[342] + v[1425] * v[344];
	v[346] = v[333] + (v[523] * v[523]);
	v[2917] = v[272] * v[346];
	v[2916] = v[270] * v[346];
	v[2919] = v[261] * v[2716];
	v[2921] = v[267] * v[2718];
	v[2298] = v[346] / v[333];
	v[1423] = v[343] / v[333];
	v[4584] = 0e0;
	v[4585] = 0e0;
	v[4586] = 0e0;
	v[4587] = 0e0;
	v[4588] = 0e0;
	v[4589] = 0e0;
	v[4590] = 0e0;
	v[4591] = 0e0;
	v[4592] = 0e0;
	v[4593] = v[2717];
	v[4594] = v[1423];
	v[4595] = 0e0;
	v[355] = v[2717] * v[336] + v[1423] * v[342] + v[2298] * v[344];
	v[2481] = v[298] - v[301] + v[309] * v[501] + v[319] * v[502] + v[329] * v[503] - v[335] * v[600] - v[345] * v[601]
		- v[355] * v[602];
	v[2477] = v[299] - v[302] + v[309] * v[504] + v[319] * v[505] + v[329] * v[506] - v[335] * v[603] - v[345] * v[604]
		- v[355] * v[605];
	v[2474] = v[300] - v[303] + v[309] * v[507] + v[319] * v[508] + v[329] * v[509] - v[335] * v[606] - v[345] * v[607]
		- v[355] * v[608];
	v[356] = GAp[0] * v[237] + GAp[1] * v[238] + GAp[2] * v[239] + v[2765] - GBp[0] * v[277] - GBp[1] * v[278] - GBp[2] * v[279];
	v[2724] = dGAp[2][0] * v[356];
	v[2719] = dGAp[2][1] * v[356];
	v[357] = GAp[0] * v[240] + GAp[1] * v[241] + GAp[2] * v[242] + v[2764] - GBp[0] * v[280] - GBp[1] * v[281] - GBp[2] * v[282];
	v[2725] = dGAp[2][0] * v[357];
	v[2720] = dGAp[2][1] * v[357];
	v[358] = GAp[0] * v[243] + GAp[1] * v[244] + GAp[2] * v[245] + v[2763] - GBp[0] * v[283] - GBp[1] * v[284] - GBp[2] * v[285];
	v[2730] = v[356] * v[459] + v[357] * v[474] + v[358] * v[489];
	v[2729] = v[356] * v[456] + v[357] * v[471] + v[358] * v[486];
	v[2728] = dGAp[2][0] * v[358];
	v[2727] = v[356] * v[460] + v[357] * v[475] + v[358] * v[490];
	v[2726] = v[356] * v[457] + v[357] * v[472] + v[358] * v[487];
	v[2723] = v[356] * v[461] + v[357] * v[476] + v[358] * v[491];
	v[2722] = v[356] * v[458] + v[357] * v[473] + v[358] * v[488];
	v[2721] = dGAp[2][1] * v[358];
	v[1264] = sqrt((v[356] * v[356]) + (v[357] * v[357]) + (v[358] * v[358]));
	v[2005] = 1e0 / Power(v[1264], 2);
	v[632] = -(v[356] * (dGBp[0][1] * v[557] + dGBp[1][1] * v[560] + dGBp[2][1] * v[563])) - v[357] * (dGBp[0][1] * v[572]
		+ dGBp[1][1] * v[575] + dGBp[2][1] * v[578]) - v[358] * (dGBp[0][1] * v[587] + dGBp[1][1] * v[590] + dGBp[2][1] * v[593]
			) + v[395] * v[602] + v[396] * v[605] + v[397] * v[608];
	v[631] = -(v[356] * (dGBp[0][1] * v[556] + dGBp[1][1] * v[559] + dGBp[2][1] * v[562])) - v[357] * (dGBp[0][1] * v[571]
		+ dGBp[1][1] * v[574] + dGBp[2][1] * v[577]) - v[358] * (dGBp[0][1] * v[586] + dGBp[1][1] * v[589] + dGBp[2][1] * v[592]
			) + v[395] * v[601] + v[396] * v[604] + v[397] * v[607];
	v[630] = -(v[356] * (dGBp[0][1] * v[555] + dGBp[1][1] * v[558] + dGBp[2][1] * v[561])) - v[357] * (dGBp[0][1] * v[570]
		+ dGBp[1][1] * v[573] + dGBp[2][1] * v[576]) - v[358] * (dGBp[0][1] * v[585] + dGBp[1][1] * v[588] + dGBp[2][1] * v[591]
			) + v[395] * v[600] + v[396] * v[603] + v[397] * v[606];
	v[626] = -(v[356] * (dGBp[0][0] * v[557] + dGBp[1][0] * v[560] + dGBp[2][0] * v[563])) - v[357] * (dGBp[0][0] * v[572]
		+ dGBp[1][0] * v[575] + dGBp[2][0] * v[578]) - v[358] * (dGBp[0][0] * v[587] + dGBp[1][0] * v[590] + dGBp[2][0] * v[593]
			) + v[391] * v[602] + v[392] * v[605] + v[393] * v[608];
	v[625] = -(v[356] * (dGBp[0][0] * v[556] + dGBp[1][0] * v[559] + dGBp[2][0] * v[562])) - v[357] * (dGBp[0][0] * v[571]
		+ dGBp[1][0] * v[574] + dGBp[2][0] * v[577]) - v[358] * (dGBp[0][0] * v[586] + dGBp[1][0] * v[589] + dGBp[2][0] * v[592]
			) + v[391] * v[601] + v[392] * v[604] + v[393] * v[607];
	v[624] = -(v[356] * (dGBp[0][0] * v[555] + dGBp[1][0] * v[558] + dGBp[2][0] * v[561])) - v[357] * (dGBp[0][0] * v[570]
		+ dGBp[1][0] * v[573] + dGBp[2][0] * v[576]) - v[358] * (dGBp[0][0] * v[585] + dGBp[1][0] * v[588] + dGBp[2][0] * v[591]
			) + v[391] * v[600] + v[392] * v[603] + v[393] * v[606];
	v[617] = dGAp[0][1] * v[2722] + dGAp[1][1] * v[2723] + v[2719] * v[464] + v[2720] * v[479] + v[2721] * v[494]
		+ v[387] * v[503] + v[388] * v[506] + v[389] * v[509];
	v[616] = dGAp[0][1] * v[2726] + dGAp[1][1] * v[2727] + v[2719] * v[463] + v[2720] * v[478] + v[2721] * v[493]
		+ v[387] * v[502] + v[388] * v[505] + v[389] * v[508];
	v[615] = dGAp[0][1] * v[2729] + dGAp[1][1] * v[2730] + v[2719] * v[462] + v[2720] * v[477] + v[2721] * v[492]
		+ v[387] * v[501] + v[388] * v[504] + v[389] * v[507];
	v[611] = dGAp[0][0] * v[2722] + dGAp[1][0] * v[2723] + v[2724] * v[464] + v[2725] * v[479] + v[2728] * v[494]
		+ v[383] * v[503] + v[384] * v[506] + v[385] * v[509];
	v[610] = dGAp[0][0] * v[2726] + dGAp[1][0] * v[2727] + v[2724] * v[463] + v[2725] * v[478] + v[2728] * v[493]
		+ v[383] * v[502] + v[384] * v[505] + v[385] * v[508];
	v[609] = dGAp[0][0] * v[2729] + dGAp[1][0] * v[2730] + v[2724] * v[462] + v[2725] * v[477] + v[2728] * v[492]
		+ v[383] * v[501] + v[384] * v[504] + v[385] * v[507];
	v[359] = GAi[0] * QAi[0][0] + GAi[1] * QAi[0][1] + GAi[2] * QAi[0][2] - GBi[0] * QBi[0][0] - GBi[1] * QBi[0][1]
		- GBi[2] * QBi[0][2] + xAi[0] - xBi[0];
	v[360] = GAi[0] * QAi[1][0] + GAi[1] * QAi[1][1] + GAi[2] * QAi[1][2] - GBi[0] * QBi[1][0] - GBi[1] * QBi[1][1]
		- GBi[2] * QBi[1][2] + xAi[1] - xBi[1];
	v[361] = GAi[0] * QAi[2][0] + GAi[1] * QAi[2][1] + GAi[2] * QAi[2][2] - GBi[0] * QBi[2][0] - GBi[1] * QBi[2][1]
		- GBi[2] * QBi[2][2] + xAi[2] - xBi[2];
	if (v[1264] > 0.1e-7) { v01 = 1e0 / v[1264]; v02 = (-(v01 / v[1264])); v03 = (2e0*v01) / Power(v[1264], 2); }
	else {
		v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[1264])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[1264])*
			(0.2399999997e10 - 0.1199999994e18*v[1264] - 0.3e17*(v[1264] * v[1264]))));
		v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[1264] + 0.6e25*Power(v[1264], 3)
			+ 0.1799999982e26*(v[1264] * v[1264]));
		v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[1264] - 0.3e17*(v[1264] * v[1264]));
	};
	v[368] = v03;
	v[369] = v02;
	v[370] = v01;
	v[371] = v[356] * v[370];
	v[2731] = (*cn)*v[371];
	v[2476] = v[2477] * v[371];
	v[2471] = v[2474] * v[371];
	v[1919] = v[2731] * v[309];
	v[1917] = v[2731] * v[329];
	v[1915] = v[2731] * v[319];
	v[1903] = -(v[2731] * v[335]);
	v[1901] = -(v[2731] * v[355]);
	v[1898] = -(v[2731] * v[345]);
	v[1880] = -(v[371] * v[602]);
	v[1879] = -(v[371] * v[601]);
	v[1878] = -(v[371] * v[600]);
	v[1877] = v[371] * v[503];
	v[1876] = v[371] * v[502];
	v[1875] = v[371] * v[501];
	v[2746] = (v[298] - v[301])*v[371];
	v[749] = (v[371] * v[371]);
	v[372] = v[357] * v[370];
	v[2483] = v[2481] * v[372];
	v[2472] = v[2474] * v[372];
	v[2288] = v[2731] * v[372];
	v[1896] = -(v[372] * v[605]);
	v[2758] = v[1880] + v[1896];
	v[1894] = -(v[372] * v[604]);
	v[2757] = v[1879] + v[1894];
	v[1892] = -(v[372] * v[603]);
	v[2756] = v[1878] + v[1892];
	v[1890] = v[372] * v[506];
	v[2755] = v[1877] + v[1890];
	v[1888] = v[372] * v[505];
	v[2754] = v[1876] + v[1888];
	v[1886] = v[372] * v[504];
	v[2753] = v[1875] + v[1886];
	v[2738] = (v[299] - v[302])*v[372];
	v[2771] = v[2738] + v[2746];
	v[2473] = v[2771] + v[2753] * v[309] + v[2754] * v[319] + v[2755] * v[329] + v[2756] * v[335] + v[2757] * v[345]
		+ v[2758] * v[355];
	v[751] = (v[372] * v[372]);
	v[373] = v[358] * v[370];
	v[3185] = v[2473] + 2e0*v[2474] * v[373];
	v[2737] = -(v[373] * v[608]);
	v[2752] = v[1880] + v[2737];
	v[2744] = -v[1896] - v[2737];
	v[2736] = -(v[373] * v[607]);
	v[2751] = v[1879] + v[2736];
	v[2743] = -v[1894] - v[2736];
	v[2735] = -(v[373] * v[606]);
	v[2750] = v[1878] + v[2735];
	v[2742] = -v[1892] - v[2735];
	v[2734] = v[373] * v[509];
	v[2749] = v[1877] + v[2734];
	v[2741] = v[1890] + v[2734];
	v[2733] = v[373] * v[508];
	v[2748] = v[1876] + v[2733];
	v[2740] = v[1888] + v[2733];
	v[2732] = v[373] * v[507];
	v[2747] = v[1875] + v[2732];
	v[2739] = v[1886] + v[2732];
	v[2484] = v[2481] * v[373];
	v[2479] = v[2477] * v[373];
	v[2290] = v[2731] * v[373];
	v[2002] = v[2747] * v[372] + v[504] * v[751];
	v[2001] = v[2739] * v[371] + v[501] * v[749];
	v[1998] = v[2748] * v[372] + v[505] * v[751];
	v[1997] = v[2740] * v[371] + v[502] * v[749];
	v[1994] = v[2749] * v[372] + v[506] * v[751];
	v[1993] = v[2741] * v[371] + v[503] * v[749];
	v[1990] = v[2750] * v[372] - v[603] * v[751];
	v[1989] = -(v[2742] * v[371]) - v[600] * v[749];
	v[1986] = v[2751] * v[372] - v[604] * v[751];
	v[1985] = -(v[2743] * v[371]) - v[601] * v[749];
	v[1982] = v[2752] * v[372] - v[605] * v[751];
	v[1981] = -(v[2744] * v[371]) - v[602] * v[749];
	v[1909] = -((*cn)*v[372] * v[373]);
	v[2745] = (v[300] - v[303])*v[373];
	v[2770] = v[2745] + v[2746];
	v[2769] = v[2738] + v[2745];
	v[2482] = v[2769] + v[2739] * v[309] + v[2740] * v[319] + v[2741] * v[329] - v[2742] * v[335] - v[2743] * v[345]
		- v[2744] * v[355];
	v[3187] = v[2482] + 2e0*v[2481] * v[371];
	v[2478] = v[2770] + v[2747] * v[309] + v[2748] * v[319] + v[2749] * v[329] + v[2750] * v[335] + v[2751] * v[345]
		+ v[2752] * v[355];
	v[3186] = v[2478] + 2e0*v[2477] * v[372];
	v[753] = (v[373] * v[373]);
	v[2003] = v[2753] * v[373] + v[507] * v[753];
	v[1999] = v[2754] * v[373] + v[508] * v[753];
	v[1995] = v[2755] * v[373] + v[509] * v[753];
	v[1991] = v[2756] * v[373] - v[606] * v[753];
	v[1987] = v[2757] * v[373] - v[607] * v[753];
	v[1983] = v[2758] * v[373] - v[608] * v[753];
	v[374] = sqrt((v[359] * v[359]) + (v[360] * v[360]) + (v[361] * v[361]));
	if (v[374] > 0.1e-7) { v04 = 1e0 / v[374]; v05 = (-(v04 / v[374])); v06 = (2e0*v04) / Power(v[374], 2); }
	else {
		v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[374])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[374])*(0.2399999997e10
			- 0.1199999994e18*v[374] - 0.3e17*(v[374] * v[374]))));
		v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[374] + 0.6e25*Power(v[374], 3)
			+ 0.1799999982e26*(v[374] * v[374]));
		v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[374] - 0.3e17*(v[374] * v[374]));
	};
	v[379] = v04;
	v[380] = v[359] * v[379];
	v[381] = v[360] * v[379];
	v[382] = v[361] * v[379];
	v[652] = -(invH[0][0] * v[609]) - invH[0][1] * v[615] - invH[0][2] * v[621] - invH[0][3] * v[627];
	v[653] = -(invH[0][0] * v[610]) - invH[0][1] * v[616] - invH[0][2] * v[622] - invH[0][3] * v[628];
	v[654] = -(invH[0][0] * v[611]) - invH[0][1] * v[617] - invH[0][2] * v[623] - invH[0][3] * v[629];
	v[661] = -(invH[0][0] * v[612]) - invH[0][1] * v[618] - invH[0][2] * v[624] - invH[0][3] * v[630];
	v[662] = -(invH[0][0] * v[613]) - invH[0][1] * v[619] - invH[0][2] * v[625] - invH[0][3] * v[631];
	v[663] = -(invH[0][0] * v[614]) - invH[0][1] * v[620] - invH[0][2] * v[626] - invH[0][3] * v[632];
	v[1162] = -(v[309] * v[652]) - v[319] * v[653] - v[329] * v[654] + v[2759] * v[655] + v[2760] * v[657] + v[2761] * v[659]
		- v[335] * v[661] - v[345] * v[662] - v[355] * v[663];
	v[3732] = v[655];
	v[3733] = v[657];
	v[3734] = v[659];
	v[3735] = v[652];
	v[3736] = v[653];
	v[3737] = v[654];
	v[3738] = -v[655];
	v[3739] = -v[657];
	v[3740] = -v[659];
	v[3741] = v[661];
	v[3742] = v[662];
	v[3743] = v[663];
	v[667] = -(invH[1][0] * v[609]) - invH[1][1] * v[615] - invH[1][2] * v[621] - invH[1][3] * v[627];
	v[668] = -(invH[1][0] * v[610]) - invH[1][1] * v[616] - invH[1][2] * v[622] - invH[1][3] * v[628];
	v[669] = -(invH[1][0] * v[611]) - invH[1][1] * v[617] - invH[1][2] * v[623] - invH[1][3] * v[629];
	v[676] = -(invH[1][0] * v[612]) - invH[1][1] * v[618] - invH[1][2] * v[624] - invH[1][3] * v[630];
	v[677] = -(invH[1][0] * v[613]) - invH[1][1] * v[619] - invH[1][2] * v[625] - invH[1][3] * v[631];
	v[678] = -(invH[1][0] * v[614]) - invH[1][1] * v[620] - invH[1][2] * v[626] - invH[1][3] * v[632];
	v[1164] = -(v[309] * v[667]) - v[319] * v[668] - v[329] * v[669] + v[2759] * v[670] + v[2760] * v[672] + v[2761] * v[674]
		- v[335] * v[676] - v[345] * v[677] - v[355] * v[678];
	v[3744] = v[670];
	v[3745] = v[672];
	v[3746] = v[674];
	v[3747] = v[667];
	v[3748] = v[668];
	v[3749] = v[669];
	v[3750] = -v[670];
	v[3751] = -v[672];
	v[3752] = -v[674];
	v[3753] = v[676];
	v[3754] = v[677];
	v[3755] = v[678];
	v[682] = -(invH[2][0] * v[609]) - invH[2][1] * v[615] - invH[2][2] * v[621] - invH[2][3] * v[627];
	v[683] = -(invH[2][0] * v[610]) - invH[2][1] * v[616] - invH[2][2] * v[622] - invH[2][3] * v[628];
	v[684] = -(invH[2][0] * v[611]) - invH[2][1] * v[617] - invH[2][2] * v[623] - invH[2][3] * v[629];
	v[691] = -(invH[2][0] * v[612]) - invH[2][1] * v[618] - invH[2][2] * v[624] - invH[2][3] * v[630];
	v[692] = -(invH[2][0] * v[613]) - invH[2][1] * v[619] - invH[2][2] * v[625] - invH[2][3] * v[631];
	v[693] = -(invH[2][0] * v[614]) - invH[2][1] * v[620] - invH[2][2] * v[626] - invH[2][3] * v[632];
	v[1158] = v[309] * v[682] + v[319] * v[683] + v[329] * v[684] - v[2759] * v[685] - v[2760] * v[687] - v[2761] * v[689]
		+ v[335] * v[691] + v[345] * v[692] + v[355] * v[693];
	v[3720] = v[685];
	v[3721] = v[687];
	v[3722] = v[689];
	v[3723] = v[682];
	v[3724] = v[683];
	v[3725] = v[684];
	v[3726] = -v[685];
	v[3727] = -v[687];
	v[3728] = -v[689];
	v[3729] = v[691];
	v[3730] = v[692];
	v[3731] = v[693];
	v[1096] = v[385] * v[655] + v[389] * v[670] - v[393] * v[685] - v[397] * v[700];
	v[1095] = v[384] * v[655] + v[388] * v[670] - v[392] * v[685] - v[396] * v[700];
	v[1094] = v[383] * v[655] + v[387] * v[670] - v[391] * v[685] - v[395] * v[700];
	v[1100] = v[385] * v[657] + v[389] * v[672] - v[393] * v[687] - v[397] * v[702];
	v[1099] = v[384] * v[657] + v[388] * v[672] - v[392] * v[687] - v[396] * v[702];
	v[1098] = v[383] * v[657] + v[387] * v[672] - v[391] * v[687] - v[395] * v[702];
	v[1104] = v[385] * v[659] + v[389] * v[674] - v[393] * v[689] - v[397] * v[704];
	v[1103] = v[384] * v[659] + v[388] * v[674] - v[392] * v[689] - v[396] * v[704];
	v[1102] = v[383] * v[659] + v[387] * v[674] - v[391] * v[689] - v[395] * v[704];
	v[697] = -(invH[3][0] * v[609]) - invH[3][1] * v[615] - invH[3][2] * v[621] - invH[3][3] * v[627];
	v[1108] = v[385] * v[652] + v[389] * v[667] - v[393] * v[682] - v[397] * v[697];
	v[1107] = v[384] * v[652] + v[388] * v[667] - v[392] * v[682] - v[396] * v[697];
	v[1106] = v[383] * v[652] + v[387] * v[667] - v[391] * v[682] - v[395] * v[697];
	v[698] = -(invH[3][0] * v[610]) - invH[3][1] * v[616] - invH[3][2] * v[622] - invH[3][3] * v[628];
	v[1112] = v[385] * v[653] + v[389] * v[668] - v[393] * v[683] - v[397] * v[698];
	v[1111] = v[384] * v[653] + v[388] * v[668] - v[392] * v[683] - v[396] * v[698];
	v[1110] = v[383] * v[653] + v[387] * v[668] - v[391] * v[683] - v[395] * v[698];
	v[699] = -(invH[3][0] * v[611]) - invH[3][1] * v[617] - invH[3][2] * v[623] - invH[3][3] * v[629];
	v[1116] = v[385] * v[654] + v[389] * v[669] - v[393] * v[684] - v[397] * v[699];
	v[1115] = v[384] * v[654] + v[388] * v[669] - v[392] * v[684] - v[396] * v[699];
	v[1114] = v[383] * v[654] + v[387] * v[669] - v[391] * v[684] - v[395] * v[699];
	v[706] = -(invH[3][0] * v[612]) - invH[3][1] * v[618] - invH[3][2] * v[624] - invH[3][3] * v[630];
	v[1132] = v[385] * v[661] + v[389] * v[676] - v[393] * v[691] - v[397] * v[706];
	v[1131] = v[384] * v[661] + v[388] * v[676] - v[392] * v[691] - v[396] * v[706];
	v[1130] = v[383] * v[661] + v[387] * v[676] - v[391] * v[691] - v[395] * v[706];
	v[707] = -(invH[3][0] * v[613]) - invH[3][1] * v[619] - invH[3][2] * v[625] - invH[3][3] * v[631];
	v[1136] = v[385] * v[662] + v[389] * v[677] - v[393] * v[692] - v[397] * v[707];
	v[1135] = v[384] * v[662] + v[388] * v[677] - v[392] * v[692] - v[396] * v[707];
	v[1134] = v[383] * v[662] + v[387] * v[677] - v[391] * v[692] - v[395] * v[707];
	v[708] = -(invH[3][0] * v[614]) - invH[3][1] * v[620] - invH[3][2] * v[626] - invH[3][3] * v[632];
	v[1160] = v[309] * v[697] + v[319] * v[698] + v[329] * v[699] - v[2759] * v[700] - v[2760] * v[702] - v[2761] * v[704]
		+ v[335] * v[706] + v[345] * v[707] + v[355] * v[708];
	v[1140] = v[385] * v[663] + v[389] * v[678] - v[393] * v[693] - v[397] * v[708];
	v[1139] = v[384] * v[663] + v[388] * v[678] - v[392] * v[693] - v[396] * v[708];
	v[1138] = v[383] * v[663] + v[387] * v[678] - v[391] * v[693] - v[395] * v[708];
	v[3708] = v[700];
	v[3709] = v[702];
	v[3710] = v[704];
	v[3711] = v[697];
	v[3712] = v[698];
	v[3713] = v[699];
	v[3714] = -v[700];
	v[3715] = -v[702];
	v[3716] = -v[704];
	v[3717] = v[706];
	v[3718] = v[707];
	v[3719] = v[708];
	b709 = sqrt(Power(v[372] * v[380] - v[371] * v[381], 2) + Power(-(v[373] * v[380]) + v[371] * v[382], 2) + Power
	(v[373] * v[381] - v[372] * v[382], 2)) > 0.1e-7;
	if (b709) {
		v[711] = v[373] * v[381] - v[372] * v[382];
		v[712] = -(v[373] * v[380]) + v[371] * v[382];
		v[713] = v[372] * v[380] - v[371] * v[381];
		v[714] = sqrt((v[711] * v[711]) + (v[712] * v[712]) + (v[713] * v[713]));
		v[1694] = 1e0 / Power(v[714], 2);
		v[1257] = v[714];
		v[1705] = 1e0 - (v[1257] * v[1257]);
		v[3065] = 1e0 / Power(v[1705], 0.15e1);
		v[1700] = 1e0 / sqrt(v[1705]);
		v[1256] = asin(v[1257]) / 2e0;
		v[3064] = tan(v[1256]);
		v[1699] = 1e0 / Power(cos(v[1256]), 2);
		v[2828] = v[1699] * v[1700];
		v[716] = 2e0*tan(v[1256]);
		if (v[714] > 0.1e-7) { v07 = 1e0 / v[714]; v08 = (-(v07 / v[714])); v09 = (2e0*v07) / Power(v[714], 2); }
		else {
			v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[714])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[714])*
				(0.2399999997e10 - 0.1199999994e18*v[714] - 0.3e17*(v[714] * v[714]))));
			v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[714] + 0.6e25*Power(v[714], 3)
				+ 0.1799999982e26*(v[714] * v[714]));
			v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[714] - 0.3e17*(v[714] * v[714]));
		};
		v[720] = v09;
		v[721] = v08;
		v[722] = v07;
		v[3063] = v[716] * v[721] + v[2828] * v[722];
		v[2762] = v[716] * v[722];
		v[723] = v[2762] * v[711];
		v[734] = (v[723] * v[723]);
		v[724] = v[2762] * v[712];
		v[732] = (v[723] * v[724]) / 2e0;
		v[727] = (v[724] * v[724]);
		v[1240] = -v[727] - v[734];
		v[725] = v[2762] * v[713];
		v[1235] = v[725] + v[732];
		v[1233] = -v[725] + v[732];
		v[739] = (v[724] * v[725]) / 2e0;
		v[1239] = v[723] + v[739];
		v[1237] = -v[723] + v[739];
		v[737] = (v[723] * v[725]) / 2e0;
		v[1238] = -v[724] + v[737];
		v[1234] = v[724] + v[737];
		v[728] = (v[725] * v[725]);
		v[1244] = 4e0 + v[727] + v[728] + v[734];
		v[3066] = 1e0 / Power(v[1244], 3);
		v[1962] = 1e0 / Power(v[1244], 2);
		v[1236] = -v[728] - v[734];
		v[1232] = -v[727] - v[728];
		v[726] = 4e0 / v[1244];
		v[729] = 1e0 + (v[1232] * v[726]) / 2e0;
		v[730] = v[1233] * v[726];
		v[731] = v[1234] * v[726];
		v[733] = v[1235] * v[726];
		v[735] = 1e0 + (v[1236] * v[726]) / 2e0;
		v[736] = v[1237] * v[726];
		v[738] = v[1238] * v[726];
		v[740] = v[1239] * v[726];
		v[741] = 1e0 + (v[1240] * v[726]) / 2e0;
	}
	else {
		v[729] = 1e0;
		v[730] = 0e0;
		v[731] = 0e0;
		v[733] = 0e0;
		v[735] = 1e0;
		v[736] = 0e0;
		v[738] = 0e0;
		v[740] = 0e0;
		v[741] = 1e0;
	};
	if ((*previouscontact)) {
		v[1198] = 1e0 - v[753];
		v[1196] = 1e0 - v[751];
		v[1194] = 1e0 - v[749];
		v[746] = GAi[0] * v[243] + GAi[1] * v[244] + GAi[2] * v[245] + v[2763] - GBi[0] * v[283] - GBi[1] * v[284] - GBi[2] * v[285]
			+ gti[0] * v[738] + gti[1] * v[740] + gti[2] * v[741];
		v[2766] = v[373] * v[746];
		v[745] = GAi[0] * v[240] + GAi[1] * v[241] + GAi[2] * v[242] + v[2764] - GBi[0] * v[280] - GBi[1] * v[281] - GBi[2] * v[282]
			+ gti[0] * v[733] + gti[1] * v[735] + gti[2] * v[736];
		v[2768] = v[372] * v[745];
		v[2789] = v[2766] + v[2768];
		v[744] = GAi[0] * v[237] + GAi[1] * v[238] + GAi[2] * v[239] + v[2765] - GBi[0] * v[277] - GBi[1] * v[278] - GBi[2] * v[279]
			+ gti[0] * v[729] + gti[1] * v[730] + gti[2] * v[731];
		v[2767] = -(v[371] * v[744]);
		v[2788] = -v[2766] + v[2767];
		v[2787] = v[2767] - v[2768];
		v[743] = -(v[2789] * v[371]) + v[1194] * v[744];
		v[747] = v[2788] * v[372] + v[1196] * v[745];
		v[748] = v[2787] * v[373] + v[1198] * v[746];
	}
	else {
		v[743] = 0e0;
		v[747] = 0e0;
		v[748] = 0e0;
	};
	v[758] = (*epsn)*v[356];
	v[759] = (*epsn)*v[357];
	v[760] = (*epsn)*v[358];
	v[761] = (*cn)*(v[2001] * v[309] + v[1997] * v[319] + v[1993] * v[329] + v[1989] * v[335] + v[1985] * v[345]
		+ v[1981] * v[355] + v[2769] * v[371] - v[2759] * v[749]);
	v[762] = (*cn)*(v[2002] * v[309] + v[1998] * v[319] + v[1994] * v[329] + v[1990] * v[335] + v[1986] * v[345]
		+ v[1982] * v[355] + v[2770] * v[372] - v[2760] * v[751]);
	v[763] = (*cn)*(v[2003] * v[309] + v[1999] * v[319] + v[1995] * v[329] + v[1991] * v[335] + v[1987] * v[345]
		+ v[1983] * v[355] + v[2771] * v[373] - v[2761] * v[753]);
	v[764] = v[758] + v[761];
	v[765] = v[759] + v[762];
	v[766] = v[760] + v[763];
	v[1816] = (v[764] * v[764]) + (v[765] * v[765]) + (v[766] * v[766]);
	v[767] = (*epst)*v[743];
	v[768] = (*epst)*v[747];
	v[769] = (*epst)*v[748];
	v[773] = -((*ct)*(v[1094] * v[298] + v[1098] * v[299] + v[1102] * v[300] - v[1094] * v[301] - v[1098] * v[302]
		- v[1102] * v[303] + v[1106] * v[309] + v[1110] * v[319] + v[1114] * v[329] + v[1130] * v[335] + v[1134] * v[345]
		+ v[1138] * v[355])) + v[767];
	v[774] = -((*ct)*(v[1095] * v[298] + v[1099] * v[299] + v[1103] * v[300] - v[1095] * v[301] - v[1099] * v[302]
		- v[1103] * v[303] + v[1107] * v[309] + v[1111] * v[319] + v[1115] * v[329] + v[1131] * v[335] + v[1135] * v[345]
		+ v[1139] * v[355])) + v[768];
	v[775] = -((*ct)*(v[1096] * v[298] + v[1100] * v[299] + v[1104] * v[300] - v[1096] * v[301] - v[1100] * v[302]
		- v[1104] * v[303] + v[1108] * v[309] + v[1112] * v[319] + v[1116] * v[329] + v[1132] * v[335] + v[1136] * v[345]
		+ v[1140] * v[355])) + v[769];
	v[1813] = (v[773] * v[773]) + (v[774] * v[774]) + (v[775] * v[775]);
	if ((*stick)) {
		b777 = sqrt((v[773] * v[773]) + (v[774] * v[774]) + (v[775] * v[775])) <= (*mus)*sqrt((v[764] * v[764]) +
			(v[765] * v[765]) + (v[766] * v[766]));
		if (b777) {
			v[779] = v[773];
			v[780] = v[774];
			v[781] = v[775];
			v[782] = 1e0;
		}
		else {
			v[1815] = 1e0 / sqrt(v[1813]);
			v[1817] = 1e0 / sqrt(v[1816]);
			v[792] = sqrt(v[1816]);
			v[783] = sqrt(v[1813]);
			if (v[783] > 0.1e-5) { v010 = 1e0 / v[783]; v011 = (-(v010 / v[783])); v012 = (2e0*v010) / Power(v[783], 2); }
			else {
				v010 = (24000000e0 - (-1e0 + 1000000e0*v[783])*(71999994e0 - 0.71999982e14*v[783] + 0.6e19*Power(v[783], 3)
					+ 0.23999982e20*(v[783] * v[783]))) / 24e0;
				v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[783] + 0.6e19*Power(v[783], 3) + 0.17999982e20*
					(v[783] * v[783]));
				v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[783] - 0.3e13*(v[783] * v[783]));
			};
			v[787] = v011;
			v[788] = v010;
			v[1814] = (*mud)*v[788] * v[792];
			v[779] = v[1814] * v[773];
			v[780] = v[1814] * v[774];
			v[781] = v[1814] * v[775];
			v[782] = 0e0;
		};
		if (sqrt((v[767] * v[767]) + (v[768] * v[768]) + (v[769] * v[769])) > (*mus)*sqrt((v[764] * v[764]) +
			(v[765] * v[765]) + (v[766] * v[766]))) {
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
			v[801] = sqrt((v[767] * v[767]) + (v[768] * v[768]) + (v[769] * v[769]));
			if (v[801] > 0.1e-5) { v016 = 1e0 / v[801]; v017 = (-(v016 / v[801])); v018 = (2e0*v016) / Power(v[801], 2); }
			else {
				v016 = (24000000e0 - (-1e0 + 1000000e0*v[801])*(71999994e0 - 0.71999982e14*v[801] + 0.6e19*Power(v[801], 3)
					+ 0.23999982e20*(v[801] * v[801]))) / 24e0;
				v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[801] + 0.6e19*Power(v[801], 3) + 0.17999982e20*
					(v[801] * v[801]));
				v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[801] - 0.3e13*(v[801] * v[801]));
			};
			v[808] = -((*mud)*v013*v016*sqrt(v[1816]));
			v[807] = v[743] + v[767] * v[808];
			v[809] = v[747] + v[768] * v[808];
			v[810] = v[748] + v[769] * v[808];
		}
		else {
			v[807] = 0e0;
			v[809] = 0e0;
			v[810] = 0e0;
		};
	}
	else {
		b811 = sqrt((v[773] * v[773]) + (v[774] * v[774]) + (v[775] * v[775])) <= (*mud)*sqrt((v[764] * v[764]) +
			(v[765] * v[765]) + (v[766] * v[766]));
		if (b811) {
			v[779] = v[773];
			v[780] = v[774];
			v[781] = v[775];
			v[782] = 1e0;
		}
		else {
			v[1823] = 1e0 / sqrt(v[1813]);
			v[1827] = 1e0 / sqrt(v[1816]);
			v[822] = sqrt(v[1816]);
			v[2854] = (*mud)*v[822];
			v[813] = sqrt(v[1813]);
			if (v[813] > 0.1e-5) { v019 = 1e0 / v[813]; v020 = (-(v019 / v[813])); v021 = (2e0*v019) / Power(v[813], 2); }
			else {
				v019 = (24000000e0 - (-1e0 + 1000000e0*v[813])*(71999994e0 - 0.71999982e14*v[813] + 0.6e19*Power(v[813], 3)
					+ 0.23999982e20*(v[813] * v[813]))) / 24e0;
				v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[813] + 0.6e19*Power(v[813], 3) + 0.17999982e20*
					(v[813] * v[813]));
				v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[813] - 0.3e13*(v[813] * v[813]));
			};
			v[817] = v020;
			v[818] = v019;
			v[2772] = (*mud)*v[818] * v[822];
			v[779] = v[2772] * v[773];
			v[780] = v[2772] * v[774];
			v[781] = v[2772] * v[775];
			v[782] = 0e0;
		};
		if (sqrt((v[767] * v[767]) + (v[768] * v[768]) + (v[769] * v[769])) > (*mud)*sqrt((v[764] * v[764]) +
			(v[765] * v[765]) + (v[766] * v[766]))) {
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
			v[831] = sqrt((v[767] * v[767]) + (v[768] * v[768]) + (v[769] * v[769]));
			if (v[831] > 0.1e-5) { v025 = 1e0 / v[831]; v026 = (-(v025 / v[831])); v027 = (2e0*v025) / Power(v[831], 2); }
			else {
				v025 = (24000000e0 - (-1e0 + 1000000e0*v[831])*(71999994e0 - 0.71999982e14*v[831] + 0.6e19*Power(v[831], 3)
					+ 0.23999982e20*(v[831] * v[831]))) / 24e0;
				v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[831] + 0.6e19*Power(v[831], 3) + 0.17999982e20*
					(v[831] * v[831]));
				v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[831] - 0.3e13*(v[831] * v[831]));
			};
			v[837] = -((*mud)*v022*v025*sqrt(v[1816]));
			v[807] = v[743] + v[767] * v[837];
			v[809] = v[747] + v[768] * v[837];
			v[810] = v[748] + v[769] * v[837];
		}
		else {
			v[807] = 0e0;
			v[809] = 0e0;
			v[810] = 0e0;
		};
	};
	fn[0] = v[764];
	fn[1] = v[765];
	fn[2] = v[766];
	ft[0] = v[779];
	ft[1] = v[780];
	ft[2] = v[781];
	(*stickupdated) = v[782];
	gtpupdated[0] = v[743] - v[807];
	gtpupdated[1] = v[747] - v[809];
	gtpupdated[2] = v[748] - v[810];
	v[849] = v[760] * v[852];
	v[850] = v[760] * v[854];
	v[902] = v[258] * v[850];
	v[851] = v[760] * v[856];
	v[908] = v[258] * v[851];
	v[853] = v[759] * v[852];
	v[903] = v[258] * v[853];
	v[855] = v[759] * v[854];
	v[905] = -(v[258] * v[855]) / 2e0;
	v[857] = v[759] * v[856];
	v[911] = v[258] * v[857];
	v[858] = v[758] * v[852];
	v[909] = v[258] * v[858];
	v[859] = v[758] * v[854];
	v[912] = v[258] * v[859];
	v[860] = v[758] * v[856];
	v[910] = -(v[258] * v[860]) / 2e0;
	v[861] = v[760] * v[864];
	v[862] = v[760] * v[866];
	v[891] = v[218] * v[862];
	v[863] = v[760] * v[868];
	v[897] = v[218] * v[863];
	v[865] = v[759] * v[864];
	v[892] = v[218] * v[865];
	v[867] = v[759] * v[866];
	v[894] = -(v[218] * v[867]) / 2e0;
	v[869] = v[759] * v[868];
	v[900] = v[218] * v[869];
	v[870] = v[758] * v[864];
	v[898] = v[218] * v[870];
	v[871] = v[758] * v[866];
	v[901] = v[218] * v[871];
	v[872] = v[758] * v[868];
	v[899] = -(v[218] * v[872]) / 2e0;
	v[873] = (*epsn)*(v[237] * v[356] + v[240] * v[357] + v[243] * v[358]);
	v[874] = (*epsn)*(v[238] * v[356] + v[241] * v[357] + v[244] * v[358]);
	v[875] = (*epsn)*(v[239] * v[356] + v[242] * v[357] + v[245] * v[358]);
	v[3028] = ddGAp[0][1][0] * v[873] + ddGAp[1][1][0] * v[874] + ddGAp[2][1][0] * v[875];
	v[3027] = ddGAp[0][0][0] * v[873] + ddGAp[1][0][0] * v[874] + ddGAp[2][0][0] * v[875];
	v[3026] = ddGAp[0][1][1] * v[873] + ddGAp[1][1][1] * v[874] + ddGAp[2][1][1] * v[875];
	v[3025] = ddGAp[0][0][1] * v[873] + ddGAp[1][0][1] * v[874] + ddGAp[2][0][1] * v[875];
	v[876] = -((*epsn)*(v[277] * v[356] + v[280] * v[357] + v[283] * v[358]));
	v[877] = -((*epsn)*(v[279] * v[356] + v[282] * v[357] + v[285] * v[358]));
	v[878] = -((*epsn)*(v[278] * v[356] + v[281] * v[357] + v[284] * v[358]));
	v[3024] = ddGBp[0][1][0] * v[876] + ddGBp[2][1][0] * v[877] + ddGBp[1][1][0] * v[878];
	v[3023] = ddGBp[0][0][0] * v[876] + ddGBp[2][0][0] * v[877] + ddGBp[1][0][0] * v[878];
	v[3022] = ddGBp[0][1][1] * v[876] + ddGBp[2][1][1] * v[877] + ddGBp[1][1][1] * v[878];
	v[3021] = ddGBp[0][0][1] * v[876] + ddGBp[2][0][1] * v[877] + ddGBp[1][0][1] * v[878];
	v[879] = v[902] + v[903];
	v[1081] = v[879] / 2e0;
	v[1077] = (v[908] + v[909]) / 2e0;
	v[881] = v[911] + v[912];
	v[1076] = v[881] / 2e0;
	v[882] = (v[550] * v[849]) / 2e0 + v[545] * v[850] + v[541] * v[851] + v[537] * v[853] + (v[532] * v[855]) / 2e0
		+ v[528] * v[857] + v[524] * v[858] + v[518] * v[859] + (v[513] * v[860]) / 2e0;
	v[1084] = v[905] + v[910] - 4e0*v[882] * v[977];
	v[1080] = v[1084] - (v[258] * v[849]) / 2e0 - v[905];
	v[1075] = v[1080] + v[905] - v[910];
	v[883] = v[891] + v[892];
	v[1071] = v[883] / 2e0;
	v[1067] = (v[897] + v[898]) / 2e0;
	v[885] = v[900] + v[901];
	v[1066] = v[885] / 2e0;
	v[886] = (v[451] * v[861]) / 2e0 + v[446] * v[862] + v[442] * v[863] + v[438] * v[865] + (v[433] * v[867]) / 2e0
		+ v[429] * v[869] + v[425] * v[870] + v[419] * v[871] + (v[414] * v[872]) / 2e0;
	v[1074] = v[894] + v[899] - 4e0*v[886] * v[962];
	v[1070] = v[1074] - (v[218] * v[861]) / 2e0 - v[894];
	v[1065] = v[1070] + v[894] - v[899];
	v[3696] = v[758];
	v[3697] = v[759];
	v[3698] = v[760];
	v[3699] = 2e0*alphaA[0] * v[1065] + alphaA[2] * v[1067] + v[400] * v[885] + v[891] - v[892];
	v[3700] = 2e0*alphaA[1] * v[1070] + (alphaA[2] * v[883]) / 2e0 + (alphaA[0] * v[885]) / 2e0 - v[897] + v[898];
	v[3701] = alphaA[0] * v[1067] + 2e0*alphaA[2] * v[1074] + v[400] * v[883] + v[900] - v[901];
	v[3702] = -v[758];
	v[3703] = -v[759];
	v[3704] = -v[760];
	v[3705] = 2e0*alphaB[0] * v[1075] + alphaB[2] * v[1077] + v[406] * v[881] + v[902] - v[903];
	v[3706] = 2e0*alphaB[1] * v[1080] + (alphaB[2] * v[879]) / 2e0 + (alphaB[0] * v[881]) / 2e0 - v[908] + v[909];
	v[3707] = alphaB[0] * v[1077] + 2e0*alphaB[2] * v[1084] + v[406] * v[879] + v[911] - v[912];
	v[887] = dGAp[0][1] * v[873] + dGAp[1][1] * v[874] + dGAp[2][1] * v[875];
	v[888] = dGAp[0][0] * v[873] + dGAp[1][0] * v[874] + dGAp[2][0] * v[875];
	v[889] = dGBp[0][0] * v[876] + dGBp[2][0] * v[877] + dGBp[1][0] * v[878];
	v[890] = dGBp[0][1] * v[876] + dGBp[2][1] * v[877] + dGBp[1][1] * v[878];
	for (i847 = 1; i847 <= 12; i847++) {
		i2778 = (i847 == 10 ? 1 : 0);
		i2777 = (i847 == 11 ? 1 : 0);
		i2776 = (i847 == 12 ? 1 : 0);
		i2775 = (i847 == 4 ? 1 : 0);
		i2774 = (i847 == 5 ? 1 : 0);
		i2773 = (i847 == 6 ? 1 : 0);
		v[932] = v[3743 + i847];
		v[931] = v[3731 + i847];
		v[927] = v[3707 + i847];
		v[926] = v[3719 + i847];
		v[918] = v[3771 + i847];
		v[964] = -4e0*v[918] * v[962];
		v[919] = v[3783 + i847];
		v[920] = v[3795 + i847];
		v[921] = v[3807 + i847];
		v[922] = v[3819 + i847];
		v[979] = -4e0*v[922] * v[977];
		v[923] = v[3831 + i847];
		v[924] = v[3843 + i847];
		v[925] = v[3855 + i847];
		v[928] = dGBp[1][0] * v[926] + dGBp[1][1] * v[927];
		v[929] = dGBp[2][0] * v[926] + dGBp[2][1] * v[927];
		v[930] = dGBp[0][0] * v[926] + dGBp[0][1] * v[927];
		v[933] = dGAp[2][0] * v[931] + dGAp[2][1] * v[932];
		v[934] = dGAp[1][0] * v[931] + dGAp[1][1] * v[932];
		v[935] = dGAp[0][0] * v[931] + dGAp[0][1] * v[932];
		v[936] = v[3867 + i847];
		v[937] = v[3915 + i847];
		v[938] = v[3963 + i847];
		v[939] = v[3975 + i847];
		v[940] = v[4023 + i847];
		v[941] = v[4071 + i847];
		v[942] = -i2773 + v[919];
		v[944] = i2773 + v[919];
		v[945] = i2774 + v[920];
		v[947] = -i2774 + v[920];
		v[948] = -i2775 + v[921];
		v[950] = i2775 + v[921];
		v[951] = -i2776 + v[923];
		v[953] = i2776 + v[923];
		v[954] = i2777 + v[924];
		v[956] = -i2777 + v[924];
		v[957] = -i2778 + v[925];
		v[959] = i2778 + v[925];
		v[961] = v[2410] * v[918] - (v[218] * v[936]) / 2e0;
		v[963] = v[218] * v[942] + v[419] * v[964];
		v[965] = v[218] * v[945] + v[425] * v[964];
		v[966] = v[218] * v[944] + v[429] * v[964];
		v[967] = (-(v[218] * v[937]) + v[433] * v[964]) / 2e0;
		v[968] = v[218] * v[948] + v[438] * v[964];
		v[969] = v[218] * v[947] + v[442] * v[964];
		v[970] = v[758] * v[961] + v[759] * v[966] + v[760] * v[969];
		v[971] = v[218] * v[950] + v[446] * v[964];
		v[972] = v[758] * v[963] + v[759] * v[967] + v[760] * v[971];
		v[973] = (-(v[218] * v[938]) + v[451] * v[964]) / 2e0;
		v[974] = v[758] * v[965] + v[759] * v[968] + v[760] * v[973];
		v[976] = v[2421] * v[922] - (v[258] * v[939]) / 2e0;
		v[978] = v[258] * v[951] + v[518] * v[979];
		v[980] = v[258] * v[954] + v[524] * v[979];
		v[981] = v[258] * v[953] + v[528] * v[979];
		v[982] = (-(v[258] * v[940]) + v[532] * v[979]) / 2e0;
		v[983] = v[258] * v[957] + v[537] * v[979];
		v[984] = v[258] * v[956] + v[541] * v[979];
		v[985] = v[758] * v[976] + v[759] * v[981] + v[760] * v[984];
		v[986] = v[258] * v[959] + v[545] * v[979];
		v[987] = v[758] * v[978] + v[759] * v[982] + v[760] * v[986];
		v[988] = (-(v[258] * v[941]) + v[550] * v[979]) / 2e0;
		v[989] = v[758] * v[980] + v[759] * v[983] + v[760] * v[988];
		v[990] = (*epsn)*(v[4083 + i847] - v[284] * v[928] - v[285] * v[929] - v[283] * v[930] + v[245] * v[933] + v[244] * v[934]
			+ v[243] * v[935] + v[868] * v[969] + v[866] * v[971] + v[864] * v[973] + v[856] * v[984] + v[854] * v[986] + v[852] * v[988]
			);
		v[991] = (*epsn)*(v[4095 + i847] - v[281] * v[928] - v[282] * v[929] - v[280] * v[930] + v[242] * v[933] + v[241] * v[934]
			+ v[240] * v[935] + v[868] * v[966] + v[866] * v[967] + v[864] * v[968] + v[856] * v[981] + v[854] * v[982] + v[852] * v[983]
			);
		v[992] = (*epsn)*(v[4107 + i847] - v[278] * v[928] - v[279] * v[929] - v[277] * v[930] + v[239] * v[933] + v[238] * v[934]
			+ v[237] * v[935] + v[868] * v[961] + v[866] * v[963] + v[864] * v[965] + v[856] * v[976] + v[854] * v[978] + v[852] * v[980]
			);
		v[993] = -(v[760] * v[929]) - GBp[2] * v[990];
		v[994] = -(v[760] * v[928]) - GBp[1] * v[990];
		v[995] = -(v[760] * v[930]) - GBp[0] * v[990];
		v[996] = v[760] * v[933] + GAp[2] * v[990];
		v[997] = v[760] * v[934] + GAp[1] * v[990];
		v[998] = v[760] * v[935] + GAp[0] * v[990];
		v[999] = -(v[759] * v[929]) - GBp[2] * v[991];
		v[1000] = -(v[759] * v[928]) - GBp[1] * v[991];
		v[1001] = -(v[759] * v[930]) - GBp[0] * v[991];
		v[1002] = v[759] * v[933] + GAp[2] * v[991];
		v[1003] = v[759] * v[934] + GAp[1] * v[991];
		v[1004] = v[759] * v[935] + GAp[0] * v[991];
		v[1005] = -(v[758] * v[929]) - GBp[2] * v[992];
		v[1006] = -(v[758] * v[928]) - GBp[1] * v[992];
		v[1007] = -(v[758] * v[930]) - GBp[0] * v[992];
		v[1008] = v[758] * v[933] + GAp[2] * v[992];
		v[1009] = v[758] * v[934] + GAp[1] * v[992];
		v[1010] = v[758] * v[935] + GAp[0] * v[992];
		v[1011] = QBi[2][2] * v[993] + QBi[2][1] * v[994] + QBi[2][0] * v[995];
		v[1012] = QBi[1][2] * v[993] + QBi[1][1] * v[994] + QBi[1][0] * v[995];
		v[1013] = QBi[0][2] * v[993] + QBi[0][1] * v[994] + QBi[0][0] * v[995];
		v[1014] = QBi[2][1] * v[1000] + QBi[2][0] * v[1001] + QBi[2][2] * v[999];
		v[1015] = QBi[1][1] * v[1000] + QBi[1][0] * v[1001] + QBi[1][2] * v[999];
		v[1016] = QBi[0][1] * v[1000] + QBi[0][0] * v[1001] + QBi[0][2] * v[999];
		v[1017] = QBi[2][2] * v[1005] + QBi[2][1] * v[1006] + QBi[2][0] * v[1007];
		v[1018] = QBi[1][2] * v[1005] + QBi[1][1] * v[1006] + QBi[1][0] * v[1007];
		v[1019] = QBi[0][2] * v[1005] + QBi[0][1] * v[1006] + QBi[0][0] * v[1007];
		v[1020] = (v[1011] * v[258] + v[849] * v[979]) / 2e0;
		v[1021] = v[1012] * v[258] + v[850] * v[979];
		v[1022] = v[1013] * v[258] + v[851] * v[979];
		v[1023] = v[1014] * v[258] + v[853] * v[979];
		v[1025] = v[1016] * v[258] + v[857] * v[979];
		v[1026] = v[1017] * v[258] + v[858] * v[979];
		v[1027] = v[1018] * v[258] + v[859] * v[979];
		v[2782] = 8e0*v[1413] * v[882] * v[922] - 4e0*((v[1019] * v[513]) / 2e0 + v[1018] * v[518] + v[1017] * v[524]
			+ v[1016] * v[528] + (v[1015] * v[532]) / 2e0 + v[1014] * v[537] + v[1013] * v[541] + v[1012] * v[545] + (v[1011] * v[550])
			/ 2e0 - (v[860] * v[939]) / 2e0 - (v[855] * v[940]) / 2e0 - (v[849] * v[941]) / 2e0 + v[859] * v[951] + v[857] * v[953]
			+ v[858] * v[954] + v[851] * v[956] + v[853] * v[957] + v[850] * v[959])*v[977];
		v[1079] = -(v[1015] * v[258]) / 2e0 + v[2782] - (v[855] * v[979]) / 2e0;
		v[1029] = (v[1019] * v[258] + v[860] * v[979]) / 2e0;
		v[1030] = QAi[2][2] * v[996] + QAi[2][1] * v[997] + QAi[2][0] * v[998];
		v[1031] = QAi[1][2] * v[996] + QAi[1][1] * v[997] + QAi[1][0] * v[998];
		v[1032] = QAi[0][2] * v[996] + QAi[0][1] * v[997] + QAi[0][0] * v[998];
		v[1033] = QAi[2][2] * v[1002] + QAi[2][1] * v[1003] + QAi[2][0] * v[1004];
		v[1034] = QAi[1][2] * v[1002] + QAi[1][1] * v[1003] + QAi[1][0] * v[1004];
		v[1035] = QAi[0][2] * v[1002] + QAi[0][1] * v[1003] + QAi[0][0] * v[1004];
		v[1036] = QAi[2][2] * v[1008] + QAi[2][1] * v[1009] + QAi[2][0] * v[1010];
		v[1037] = QAi[1][2] * v[1008] + QAi[1][1] * v[1009] + QAi[1][0] * v[1010];
		v[1038] = QAi[0][2] * v[1008] + QAi[0][1] * v[1009] + QAi[0][0] * v[1010];
		v[1039] = (v[1030] * v[218] + v[861] * v[964]) / 2e0;
		v[1040] = v[1031] * v[218] + v[862] * v[964];
		v[1041] = v[1032] * v[218] + v[863] * v[964];
		v[1042] = v[1033] * v[218] + v[865] * v[964];
		v[1044] = v[1035] * v[218] + v[869] * v[964];
		v[1045] = v[1036] * v[218] + v[870] * v[964];
		v[1046] = v[1037] * v[218] + v[871] * v[964];
		v[2780] = 8e0*v[1411] * v[886] * v[918] - 4e0*((v[1038] * v[414]) / 2e0 + v[1037] * v[419] + v[1036] * v[425]
			+ v[1035] * v[429] + (v[1034] * v[433]) / 2e0 + v[1033] * v[438] + v[1032] * v[442] + v[1031] * v[446] + (v[1030] * v[451])
			/ 2e0 - (v[872] * v[936]) / 2e0 - (v[867] * v[937]) / 2e0 - (v[861] * v[938]) / 2e0 + v[871] * v[942] + v[869] * v[944]
			+ v[870] * v[945] + v[863] * v[947] + v[865] * v[948] + v[862] * v[950])*v[962];
		v[1069] = -(v[1034] * v[218]) / 2e0 + v[2780] - (v[867] * v[964]) / 2e0;
		v[1048] = (v[1038] * v[218] + v[872] * v[964]) / 2e0;
		v[2781] = (v[1022] + v[1026]) / 2e0;
		v[1050] = v[1021] + v[1023];
		v[1051] = v[1025] + v[1027];
		v[1052] = -(QBi[0][2] * v[985]) - QBi[1][2] * v[987] - QBi[2][2] * v[989] - v[285] * v[990] - v[282] * v[991]
			- v[279] * v[992];
		v[1053] = -(QBi[0][1] * v[985]) - QBi[1][1] * v[987] - QBi[2][1] * v[989] - v[284] * v[990] - v[281] * v[991]
			- v[278] * v[992];
		v[1054] = -(QBi[0][0] * v[985]) - QBi[1][0] * v[987] - QBi[2][0] * v[989] - v[283] * v[990] - v[280] * v[991]
			- v[277] * v[992];
		v[1055] = dGBp[2][1] * v[1052] + dGBp[1][1] * v[1053] + dGBp[0][1] * v[1054] + v[3021] * v[926] + v[3022] * v[927];
		v[1056] = dGBp[2][0] * v[1052] + dGBp[1][0] * v[1053] + dGBp[0][0] * v[1054] + v[3023] * v[926] + v[3024] * v[927];
		v[2779] = (v[1041] + v[1045]) / 2e0;
		v[1058] = v[1040] + v[1042];
		v[1059] = v[1044] + v[1046];
		v[4120] = 0e0;
		v[4121] = 0e0;
		v[4122] = 0e0;
		v[4123] = 2e0*v[1065];
		v[4124] = v[1066];
		v[4125] = v[1067];
		v[4126] = 0e0;
		v[4127] = 0e0;
		v[4128] = 0e0;
		v[4129] = 0e0;
		v[4130] = 0e0;
		v[4131] = 0e0;
		v[4132] = 0e0;
		v[4133] = 0e0;
		v[4134] = 0e0;
		v[4135] = v[1066];
		v[4136] = 2e0*v[1070];
		v[4137] = v[1071];
		v[4138] = 0e0;
		v[4139] = 0e0;
		v[4140] = 0e0;
		v[4141] = 0e0;
		v[4142] = 0e0;
		v[4143] = 0e0;
		v[4144] = 0e0;
		v[4145] = 0e0;
		v[4146] = 0e0;
		v[4147] = v[1067];
		v[4148] = v[1071];
		v[4149] = 2e0*v[1074];
		v[4150] = 0e0;
		v[4151] = 0e0;
		v[4152] = 0e0;
		v[4153] = 0e0;
		v[4154] = 0e0;
		v[4155] = 0e0;
		v[4156] = 0e0;
		v[4157] = 0e0;
		v[4158] = 0e0;
		v[4159] = 0e0;
		v[4160] = 0e0;
		v[4161] = 0e0;
		v[4162] = 0e0;
		v[4163] = 0e0;
		v[4164] = 0e0;
		v[4165] = 2e0*v[1075];
		v[4166] = v[1076];
		v[4167] = v[1077];
		v[4168] = 0e0;
		v[4169] = 0e0;
		v[4170] = 0e0;
		v[4171] = 0e0;
		v[4172] = 0e0;
		v[4173] = 0e0;
		v[4174] = 0e0;
		v[4175] = 0e0;
		v[4176] = 0e0;
		v[4177] = v[1076];
		v[4178] = 2e0*v[1080];
		v[4179] = v[1081];
		v[4180] = 0e0;
		v[4181] = 0e0;
		v[4182] = 0e0;
		v[4183] = 0e0;
		v[4184] = 0e0;
		v[4185] = 0e0;
		v[4186] = 0e0;
		v[4187] = 0e0;
		v[4188] = 0e0;
		v[4189] = v[1077];
		v[4190] = v[1081];
		v[4191] = 2e0*v[1084];
		v[4192] = v[992];
		v[4193] = v[991];
		v[4194] = v[990];
		v[4195] = v[1040] - v[1042] + 2e0*alphaA[0] * (-v[1039] + v[1069]) + alphaA[2] * v[2779] + v[1059] * v[400] + v[4119
			+ i847];
		v[4196] = -v[1041] + v[1045] + (alphaA[2] * v[1058]) / 2e0 + (alphaA[0] * v[1059]) / 2e0 + 2e0*alphaA[1] * (-v[1039]
			- v[1048] + v[2780]) + v[4131 + i847];
		v[4197] = v[1044] - v[1046] + 2e0*alphaA[2] * (-v[1048] + v[1069]) + alphaA[0] * v[2779] + v[1058] * v[400] + v[4143
			+ i847];
		v[4198] = -v[992];
		v[4199] = -v[991];
		v[4200] = -v[990];
		v[4201] = v[1021] - v[1023] + 2e0*alphaB[0] * (-v[1020] + v[1079]) + alphaB[2] * v[2781] + v[1051] * v[406] + v[4155
			+ i847];
		v[4202] = -v[1022] + v[1026] + (alphaB[2] * v[1050]) / 2e0 + (alphaB[0] * v[1051]) / 2e0 + 2e0*alphaB[1] * (-v[1020]
			- v[1029] + v[2782]) + v[4167 + i847];
		v[4203] = v[1025] - v[1027] + 2e0*alphaB[2] * (-v[1029] + v[1079]) + alphaB[0] * v[2781] + v[1050] * v[406] + v[4179
			+ i847];
		v[1060] = QAi[0][0] * v[970] + QAi[1][0] * v[972] + QAi[2][0] * v[974] + v[243] * v[990] + v[240] * v[991]
			+ v[237] * v[992];
		v[1061] = QAi[0][1] * v[970] + QAi[1][1] * v[972] + QAi[2][1] * v[974] + v[244] * v[990] + v[241] * v[991]
			+ v[238] * v[992];
		v[1062] = QAi[0][2] * v[970] + QAi[1][2] * v[972] + QAi[2][2] * v[974] + v[245] * v[990] + v[242] * v[991]
			+ v[239] * v[992];
		v[1063] = dGAp[0][1] * v[1060] + dGAp[1][1] * v[1061] + dGAp[2][1] * v[1062] + v[3025] * v[931] + v[3026] * v[932];
		v[1064] = dGAp[0][0] * v[1060] + dGAp[1][0] * v[1061] + dGAp[2][0] * v[1062] + v[3027] * v[931] + v[3028] * v[932];
		Rc[i847 - 1] += v[3695 + i847] + v[889] * v[926] + v[890] * v[927] + v[888] * v[931] + v[887] * v[932];
		for (i915 = i847; i915 <= 12; i915++) {
			v[1089] = v[1055] * v[3707 + i915] + v[1056] * v[3719 + i915] + v[1064] * v[3731 + i915] + v[1063] * v[3743 + i915]
				+ v[4191 + i915];
			Kc[i847 - 1][i915 - 1] += v[1089];
			if (i847 != i915) {
				Kc[i915 - 1][i847 - 1] += v[1089];
			}
			else {
			};
		};/* end for */
	};/* end for */
	v[1097] = -(v[1094] * v[779]) - v[1095] * v[780] - v[1096] * v[781];
	v[1101] = -(v[1098] * v[779]) - v[1099] * v[780] - v[1100] * v[781];
	v[1105] = -(v[1102] * v[779]) - v[1103] * v[780] - v[1104] * v[781];
	v[1145] = 0e0;
	v[1146] = 0e0;
	v[1147] = 0e0;
	b1148 = (*stick);
	if (b1148) {
		b1149 = b777;
		if (b1149) {
			v[1147] = 0e0;
			v[1146] = 0e0;
			v[1145] = 0e0;
		}
		else {
		};
	}
	else {
		b1150 = b811;
		if (b1150) {
			v[1147] = 0e0;
			v[1146] = 0e0;
			v[1145] = 0e0;
		}
		else {
		};
	};
	v[2786] = (*ct)*v[1145];
	v[2862] = v[1145] * v[2783] - v[779];
	v[2785] = (*ct)*v[1146];
	v[2860] = v[1146] * v[2783] - v[780];
	v[2944] = (*ct)*(-(v[1102] * v[1145]) - v[1103] * v[1146] - v[1104] * v[1147]);
	v[2943] = (*ct)*(-(v[1098] * v[1145]) - v[1099] * v[1146] - v[1100] * v[1147]);
	v[2942] = (*ct)*(-(v[1094] * v[1145]) - v[1095] * v[1146] - v[1096] * v[1147]);
	v[2784] = (*ct)*v[1147];
	v[2858] = v[1147] * v[2783] - v[781];
	v[1154] = v[1158] * v[2784];
	v[1155] = v[1160] * v[2784];
	v[1156] = v[1162] * v[2784];
	v[1157] = v[1164] * v[2784];
	v[1159] = v[1158] * v[2785];
	v[1161] = v[1160] * v[2785];
	v[1163] = v[1162] * v[2785];
	v[1165] = v[1164] * v[2785];
	v[1166] = -((*ct)*(v[1138] * v[1145] + v[1139] * v[1146] + v[1140] * v[1147]));
	v[1350] = v[1166] * v[1426];
	v[1330] = v[1166] * v[2298];
	v[1167] = -((*ct)*(v[1134] * v[1145] + v[1135] * v[1146] + v[1136] * v[1147]));
	v[1337] = v[1167] * v[2297];
	v[2794] = v[1337] + v[1166] * v[1425];
	v[2793] = v[1330] + v[1167] * v[1423];
	v[1168] = -((*ct)*(v[1130] * v[1145] + v[1131] * v[1146] + v[1132] * v[1147]));
	v[1348] = v[1168] * v[2296];
	v[2989] = v[1348] + v[1350];
	v[2795] = v[1348] + v[1167] * v[1424];
	v[2980] = v[1350] + v[2795];
	v[1339] = v[1168] * v[2715];
	v[2986] = v[1339] + v[2794];
	v[2981] = v[1337] + v[1339];
	v[1333] = v[1168] * v[2717];
	v[2990] = v[1333] + v[2793];
	v[2982] = v[1330] + v[1333];
	v[1169] = -((*ct)*(v[1114] * v[1145] + v[1115] * v[1146] + v[1116] * v[1147]));
	v[1388] = v[1169] * v[1420];
	v[1368] = v[1169] * v[2295];
	v[1170] = -((*ct)*(v[1110] * v[1145] + v[1111] * v[1146] + v[1112] * v[1147]));
	v[1375] = v[1170] * v[2294];
	v[2807] = v[1375] + v[1169] * v[1419];
	v[2806] = v[1368] + v[1170] * v[1417];
	v[1171] = -((*ct)*(v[1106] * v[1145] + v[1107] * v[1146] + v[1108] * v[1147]));
	v[1386] = v[1171] * v[2293];
	v[4220] = v[2942];
	v[4221] = v[2943];
	v[4222] = v[2944];
	v[4223] = v[1386] + v[1170] * v[2710] + v[1169] * v[2712];
	v[4224] = v[1375] + v[1169] * v[1417] + v[1171] * v[1418];
	v[4225] = v[1368] + v[1170] * v[1419] + v[1171] * v[1420];
	v[4226] = -v[2942];
	v[4227] = -v[2943];
	v[4228] = -v[2944];
	v[4229] = v[1348] + v[1167] * v[2715] + v[1166] * v[2717];
	v[4230] = v[1337] + v[1166] * v[1423] + v[1168] * v[1424];
	v[4231] = v[1330] + v[1167] * v[1425] + v[1168] * v[1426];
	v[2977] = v[1386] + v[1388];
	v[2808] = v[1386] + v[1170] * v[1418];
	v[2968] = v[1388] + v[2808];
	v[1377] = v[1171] * v[2710];
	v[2974] = v[1377] + v[2807];
	v[2969] = v[1375] + v[1377];
	v[1371] = v[1171] * v[2712];
	v[2978] = v[1371] + v[2806];
	v[2970] = v[1368] + v[1371];
	v[1172] = v[1158] * v[2786];
	v[1173] = v[1160] * v[2786];
	v[1174] = v[1162] * v[2786];
	v[1175] = v[1164] * v[2786];
	v[1176] = (*epst)*v[1147];
	v[1752] = -(v[1176] * v[373]);
	v[1177] = (*epst)*v[1146];
	v[1754] = -(v[1177] * v[372]);
	v[1751] = v[1752] + v[1754];
	v[1178] = (*epst)*v[1145];
	v[1755] = -(v[1178] * v[371]);
	v[1756] = v[1754] + v[1755];
	v[1753] = v[1752] + v[1755];
	v[1179] = 0e0;
	v[1180] = 0e0;
	v[1181] = 0e0;
	v[1182] = 0e0;
	v[1183] = 0e0;
	v[1184] = 0e0;
	v[1185] = 0e0;
	v[1186] = 0e0;
	v[1187] = 0e0;
	v[1188] = 0e0;
	v[1189] = 0e0;
	v[1190] = 0e0;
	b1191 = (*previouscontact);
	if (b1191) {
		v[1192] = -(v[1176] * v[746]);
		v[1193] = -(v[1177] * v[745]);
		v[1195] = v[1178] * v[1194] + v[1751] * v[371];
		v[1197] = v[1177] * v[1196] + v[1753] * v[372];
		v[1199] = v[1176] * v[1198] + v[1756] * v[373];
		v[1181] = v[1176] * v[2787] + v[1756] * v[746];
		v[1180] = v[1177] * v[2788] + v[1753] * v[745];
		v[1201] = -(v[1178] * v[744]);
		v[1179] = -(v[1178] * v[2789]) + v[1751] * v[744];
		v[1182] = gti[0] * v[1195];
		v[1183] = gti[1] * v[1195];
		v[1184] = gti[2] * v[1195];
		v[1204] = -v[1195];
		v[1205] = -(GBi[2] * v[1195]);
		v[1206] = -(GBi[1] * v[1195]);
		v[1207] = -(GBi[0] * v[1195]);
		v[1208] = v[1195];
		v[1209] = GAi[2] * v[1195];
		v[1210] = GAi[1] * v[1195];
		v[1211] = GAi[0] * v[1195];
		v[1185] = gti[0] * v[1197];
		v[1186] = gti[1] * v[1197];
		v[1187] = gti[2] * v[1197];
		v[1212] = -v[1197];
		v[1213] = -(GBi[2] * v[1197]);
		v[1214] = -(GBi[1] * v[1197]);
		v[1215] = -(GBi[0] * v[1197]);
		v[1216] = v[1197];
		v[1217] = GAi[2] * v[1197];
		v[1218] = GAi[1] * v[1197];
		v[1219] = GAi[0] * v[1197];
		v[1188] = gti[0] * v[1199];
		v[1189] = gti[1] * v[1199];
		v[1190] = gti[2] * v[1199];
		v[1220] = -v[1199];
		v[1221] = -(GBi[2] * v[1199]);
		v[1222] = -(GBi[1] * v[1199]);
		v[1223] = -(GBi[0] * v[1199]);
		v[1224] = v[1199];
		v[1225] = GAi[2] * v[1199];
		v[1226] = GAi[1] * v[1199];
		v[1227] = GAi[0] * v[1199];
	}
	else {
		v[1211] = 0e0;
		v[1210] = 0e0;
		v[1209] = 0e0;
		v[1219] = 0e0;
		v[1218] = 0e0;
		v[1217] = 0e0;
		v[1227] = 0e0;
		v[1226] = 0e0;
		v[1225] = 0e0;
		v[1208] = 0e0;
		v[1216] = 0e0;
		v[1224] = 0e0;
		v[1207] = 0e0;
		v[1206] = 0e0;
		v[1205] = 0e0;
		v[1215] = 0e0;
		v[1214] = 0e0;
		v[1213] = 0e0;
		v[1223] = 0e0;
		v[1222] = 0e0;
		v[1221] = 0e0;
		v[1204] = 0e0;
		v[1212] = 0e0;
		v[1220] = 0e0;
		v[1201] = 0e0;
		v[1193] = 0e0;
		v[1192] = 0e0;
	};
	b1228 = b709;
	if (b1228) {
		v[1254] = -(v[1190] * v[726]) / 2e0;
		v[1253] = -(v[1186] * v[726]) / 2e0;
		v[1252] = v[1189] * v[726];
		v[1251] = v[1187] * v[726];
		v[1247] = v[1188] * v[726];
		v[1246] = v[1184] * v[726];
		v[1243] = v[1185] * v[726];
		v[1242] = v[1183] * v[726];
		v[1229] = v[1251] + v[1252];
		v[1230] = v[1246] + v[1247];
		v[1231] = v[1242] + v[1243];
		v[1241] = (v[1182] * v[1232]) / 2e0 + v[1183] * v[1233] + v[1184] * v[1234] + v[1185] * v[1235] + (v[1186] * v[1236])
			/ 2e0 + v[1187] * v[1237] + v[1188] * v[1238] + v[1189] * v[1239] + (v[1190] * v[1240]) / 2e0;
		v[1718] = (-4e0*v[1241]) / Power(v[1244], 2) + v[1253] + v[1254];
		v[1717] = -v[1253] + v[1718] - (v[1182] * v[726]) / 2e0;
		v[1716] = v[1253] - v[1254] + v[1717];
		v[1245] = (-2e0*v[1242] + 2e0*v[1243] + v[1230] * v[723] + v[1229] * v[724] + 4e0*v[1716] * v[725]) / 2e0;
		v[1250] = (2e0*v[1246] - 2e0*v[1247] + v[1231] * v[723] + 4e0*v[1717] * v[724] + v[1229] * v[725]) / 2e0;
		v[1255] = (-2e0*v[1251] + 2e0*v[1252] + 4e0*v[1718] * v[723] + v[1231] * v[724] + v[1230] * v[725]) / 2e0;
		v[2790] = v[1255] * v[711] + v[1250] * v[712] + v[1245] * v[713];
		v[1704] = v[2790] * v[722];
		v[1701] = v[2790] * v[716];
		v[1258] = v[1701] * v[721] + v[1704] / (Power(cos(v[1256]), 2)*sqrt(v[1705]));
		v[2832] = v[1258] / v[714];
		v[2791] = v[1258] / v[714];
		v[1259] = v[1245] * v[2762] + v[2791] * v[713];
		v[1261] = v[1250] * v[2762] + v[2791] * v[712];
		v[1262] = v[1255] * v[2762] + v[2791] * v[711];
		v[1179] = v[1179] - v[1259] * v[381] + v[1261] * v[382];
		v[1181] = v[1181] - v[1261] * v[380] + v[1262] * v[381];
		v[1180] = v[1180] + v[1259] * v[380] - v[1262] * v[382];
	}
	else {
	};
	v[1181] = v[1181] + 2e0*v[1192] * v[373];
	v[1180] = v[1180] + 2e0*v[1193] * v[372];
	v[1179] = v[1179] + 2e0*v[1201] * v[371];
	v[2006] = v[1179] * v[356] + v[1180] * v[357] + v[1181] * v[358];
	v[1263] = v[2006] * v[369];
	v[2792] = v[1263] / v[1264];
	v[1265] = v[2792] * v[358] + v[1181] * v[370];
	v[1267] = v[2792] * v[357] + v[1180] * v[370];
	v[1268] = v[2792] * v[356] + v[1179] * v[370];
	v[1220] = v[1220] - v[1265];
	v[1224] = v[1224] + v[1265];
	v[1212] = v[1212] - v[1267];
	v[1216] = v[1216] + v[1267];
	v[1204] = v[1204] - v[1268];
	v[1208] = v[1208] + v[1268];
	v[1275] = v[1166] * v[1785];
	v[1276] = v[1166] * v[1784];
	v[1277] = v[1166] * v[1783];
	v[1280] = v[1167] * v[1780];
	v[1281] = v[1166] * v[342] + v[1167] * v[344];
	v[2796] = v[1281] * v[263];
	v[1284] = v[1167] * v[1778];
	v[1285] = v[1275] + v[1280];
	v[1286] = -v[1275] + v[1280];
	v[2953] = 8e0*v[1286];
	v[1288] = v[1167] * v[1779] + v[2796] / v[333];
	v[1291] = v[1168] * v[1773];
	v[1292] = v[1166] * v[336] + v[1168] * v[344];
	v[2797] = v[1292] * v[268];
	v[1293] = v[1167] * v[336] + v[1168] * v[342];
	v[2798] = v[1293] * v[272];
	v[3172] = -(v[2714] * v[2796]) - v[2798] * v[331] - v[2797] * v[332] - v[1167] * (v[2918] * v[336] + v[2922] * v[336]
		+ v[337] * v[342] + v[2923] * v[344] + v[2925] * v[344]) - v[1168] * (v[2920] * v[342] + v[2924] * v[344] + v[334] * (v[336]
			+ v[262] * v[342] + v[263] * v[344])) - v[1166] * (v[2916] * v[336] + v[2919] * v[336] + v[2917] * v[342] + v[2921] * v[342]
				+ v[344] * v[346]);
	v[1295] = v[1168] * v[1775] + v[2797] / v[333];
	v[1296] = v[1288] + v[1295];
	v[1297] = -v[1288] + v[1295];
	v[2951] = 8e0*v[1297];
	v[1299] = v[1168] * v[1774] + v[2798] / v[333];
	v[1300] = v[1276] + v[1299];
	v[1301] = -v[1276] + v[1299];
	v[3168] = v[1286] * v[405] - v[1301] * v[408] + v[1297] * v[410];
	v[2952] = -8e0*v[1301];
	v[4980] = 0e0;
	v[4981] = 0e0;
	v[4982] = 0e0;
	v[4983] = 0e0;
	v[4984] = 0e0;
	v[4985] = 0e0;
	v[4986] = 0e0;
	v[4987] = 0e0;
	v[4988] = 0e0;
	v[4989] = v[2953];
	v[4990] = v[2952];
	v[4991] = v[2951];
	v[1302] = v[1169] * v[1770];
	v[1303] = v[1169] * v[1769];
	v[1304] = v[1169] * v[1768];
	v[1307] = v[1170] * v[1765];
	v[1308] = v[1169] * v[316] + v[1170] * v[318];
	v[2809] = v[1308] * v[223];
	v[1311] = v[1170] * v[1763];
	v[1312] = v[1302] + v[1307];
	v[1313] = -v[1302] + v[1307];
	v[2947] = 8e0*v[1313];
	v[1315] = v[1170] * v[1764] + v[2809] / v[307];
	v[1318] = v[1171] * v[1758];
	v[1319] = v[1169] * v[310] + v[1171] * v[318];
	v[2810] = v[1319] * v[228];
	v[1320] = v[1170] * v[310] + v[1171] * v[316];
	v[2811] = v[1320] * v[232];
	v[3178] = -(v[2709] * v[2809]) - v[2811] * v[305] - v[2810] * v[306] - v[1170] * (v[2931] * v[310] + v[2935] * v[310]
		+ v[311] * v[316] + v[2936] * v[318] + v[2938] * v[318]) - v[1171] * (v[2933] * v[316] + v[2937] * v[318] + v[308] * (v[310]
			+ v[222] * v[316] + v[223] * v[318])) - v[1169] * (v[2929] * v[310] + v[2932] * v[310] + v[2930] * v[316] + v[2934] * v[316]
				+ v[318] * v[320]);
	v[1322] = v[1171] * v[1760] + v[2810] / v[307];
	v[1323] = v[1315] + v[1322];
	v[1324] = -v[1315] + v[1322];
	v[2945] = 8e0*v[1324];
	v[1326] = v[1171] * v[1759] + v[2811] / v[307];
	v[1327] = v[1303] + v[1326];
	v[1328] = -v[1303] + v[1326];
	v[3174] = v[1313] * v[399] - v[1328] * v[402] + v[1324] * v[404];
	v[2946] = -8e0*v[1328];
	v[4992] = 0e0;
	v[4993] = 0e0;
	v[4994] = 0e0;
	v[4995] = v[2947];
	v[4996] = v[2946];
	v[4997] = v[2945];
	v[4998] = 0e0;
	v[4999] = 0e0;
	v[5000] = 0e0;
	v[5001] = 0e0;
	v[5002] = 0e0;
	v[5003] = 0e0;
	v[1221] = dGBp[2][0] * v[1154] + dGBp[2][1] * v[1155] + v[1221] - GBp[2] * v[1265];
	v[1222] = dGBp[1][0] * v[1154] + dGBp[1][1] * v[1155] + v[1222] - GBp[1] * v[1265];
	v[1223] = dGBp[0][0] * v[1154] + dGBp[0][1] * v[1155] + v[1223] - GBp[0] * v[1265];
	v[1329] = QBi[2][2] * v[1221] + QBi[2][1] * v[1222] + QBi[2][0] * v[1223] + v[2990] * v[344];
	v[1332] = QBi[1][2] * v[1221] + QBi[1][1] * v[1222] + QBi[1][0] * v[1223] + v[1293] * v[2717] + v[2793] * v[342];
	v[1334] = QBi[0][2] * v[1221] + QBi[0][1] * v[1222] + QBi[0][0] * v[1223] + v[2982] * v[336];
	v[1213] = dGBp[2][0] * v[1159] + dGBp[2][1] * v[1161] + v[1213] - GBp[2] * v[1267];
	v[1214] = dGBp[1][0] * v[1159] + dGBp[1][1] * v[1161] + v[1214] - GBp[1] * v[1267];
	v[1215] = dGBp[0][0] * v[1159] + dGBp[0][1] * v[1161] + v[1215] - GBp[0] * v[1267];
	v[1335] = QBi[2][2] * v[1213] + QBi[2][1] * v[1214] + QBi[2][0] * v[1215] + v[1292] * v[2715] + v[2794] * v[344];
	v[1338] = QBi[1][2] * v[1213] + QBi[1][1] * v[1214] + QBi[1][0] * v[1215] + v[2986] * v[342];
	v[1340] = QBi[0][2] * v[1213] + QBi[0][1] * v[1214] + QBi[0][0] * v[1215] + v[2981] * v[336];
	v[1205] = dGBp[2][0] * v[1172] + dGBp[2][1] * v[1173] + v[1205] - GBp[2] * v[1268];
	v[1206] = dGBp[1][0] * v[1172] + dGBp[1][1] * v[1173] + v[1206] - GBp[1] * v[1268];
	v[1207] = dGBp[0][0] * v[1172] + dGBp[0][1] * v[1173] + v[1207] - GBp[0] * v[1268];
	v[1347] = QBi[2][2] * v[1205] + QBi[2][1] * v[1206] + QBi[2][0] * v[1207] + v[1281] * v[1424] + v[2989] * v[344];
	v[1349] = QBi[1][2] * v[1205] + QBi[1][1] * v[1206] + QBi[1][0] * v[1207] + v[2795] * v[342];
	v[1352] = QBi[0][2] * v[1205] + QBi[0][1] * v[1206] + QBi[0][0] * v[1207] + v[2980] * v[336];
	v[1353] = -(v[1296] * v[517]) + 2e0*v[1291] * v[520] + v[1300] * v[523];
	v[2956] = -v[1353] / 2e0;
	v[1354] = 2e0*v[1284] * v[517] - v[1296] * v[520] - v[1285] * v[523];
	v[2955] = v[1354] / 2e0;
	v[1355] = -(v[1285] * v[517]) + v[1300] * v[520] + 2e0*v[1277] * v[523];
	v[2954] = -v[1355] / 2e0;
	v[4956] = 0e0;
	v[4957] = 0e0;
	v[4958] = 0e0;
	v[4959] = 0e0;
	v[4960] = 0e0;
	v[4961] = 0e0;
	v[4962] = 0e0;
	v[4963] = 0e0;
	v[4964] = 0e0;
	v[4965] = v[2956];
	v[4966] = v[2955];
	v[4967] = v[2954];
	v[1356] = (v[1329] * v[258]) / 2e0;
	v[1357] = v[1332] * v[258];
	v[1358] = v[1334] * v[258];
	v[1359] = v[1335] * v[258];
	v[1360] = (v[1338] * v[258]) / 2e0;
	v[1361] = v[1340] * v[258];
	v[1362] = v[1347] * v[258];
	v[1363] = v[1349] * v[258];
	v[1365] = 1e0 / Power(v[333], 2);
	v[2805] = -(v[1365] * v[342]);
	v[2804] = -(v[1365] * v[344]);
	v[2803] = -(v[1365] * v[336]);
	v[2802] = -(v[1365] * v[331]);
	v[2801] = -(v[1166] * v[1365]);
	v[2800] = -(v[1365] * v[332]);
	v[2799] = -(v[1365] * v[2714]);
	v[2234] = -(v[1365] * (v[2924] + v[263] * v[334]));
	v[2233] = -(v[1365] * (v[2920] + v[262] * v[334]));
	v[2232] = -(v[1365] * (v[2918] + v[2922]));
	v[2231] = -(v[1365] * (v[2923] + v[2925]));
	v[2230] = -(v[1365] * (v[2916] + v[2919]));
	v[2229] = -(v[1365] * (v[2917] + v[2921]));
	v[2228] = -(v[1365] * v[334]);
	v[2928] = v[2228] * v[336];
	v[2227] = -(v[1365] * v[337]);
	v[2927] = v[2227] * v[342];
	v[2926] = -(v[1365] * v[344] * v[346]);
	v[2225] = -(v[1365] * v[2796]);
	v[2224] = -(v[1365] * v[2797]);
	v[2223] = -(v[1365] * v[2798]);
	v[2222] = v[1293] * v[2802];
	v[2221] = v[1292] * v[2800];
	v[2220] = v[1281] * v[2799];
	v[2099] = v[1167] * v[2799];
	v[2098] = v[2716] * v[2801];
	v[2096] = v[1168] * v[2228];
	v[3157] = v[2220] + (v[2096] + v[2098])*v[344];
	v[2912] = v[2096] + v[2099];
	v[3158] = v[2098] + v[2912];
	v[2092] = v[1168] * v[2800];
	v[2089] = v[1167] * v[2227];
	v[3156] = v[2089] + v[2092];
	v[2088] = v[2718] * v[2801];
	v[2911] = v[2088] + v[2089];
	v[3155] = v[2092] + v[2911];
	v[3154] = v[2221] + v[2911] * v[344];
	v[2084] = v[1168] * v[2802];
	v[2082] = -(v[1167] * v[1365] * v[343]);
	v[2081] = v[2801] * v[346];
	v[3153] = v[2081] + v[2084];
	v[2910] = v[2081] + v[2082];
	v[3152] = v[2222] + v[2910] * v[342];
	v[3151] = v[2084] + v[2910];
	v[2041] = v[270] * v[2803];
	v[2038] = v[265] * v[2803];
	v[2035] = v[263] * v[2804];
	v[2034] = v[262] * v[2805];
	v[2029] = v[268] * v[2804];
	v[2028] = v[267] * v[2805];
	v[3143] = v[2028] + v[2038];
	v[2899] = v[2028] + v[2029];
	v[3141] = v[2038] + v[2899];
	v[2026] = v[261] * v[2803];
	v[3139] = v[2026] + v[2035];
	v[2900] = v[2026] + v[2034];
	v[3142] = v[2035] + v[2900];
	v[2022] = v[273] * v[2804];
	v[3144] = v[2022] + v[2041];
	v[2021] = v[272] * v[2805];
	v[2897] = v[2021] + v[2022];
	v[3140] = v[2041] + v[2897];
	v[1782] = v[2926] + v[2230] * v[336] + v[2229] * v[342];
	v[1777] = v[2927] + v[2232] * v[336] + v[2231] * v[344];
	v[1772] = v[2928] + v[2233] * v[342] + v[2234] * v[344];
	v[1613] = v[263] * v[2799];
	v[1611] = v[268] * v[2800];
	v[1607] = v[272] * v[2802];
	v[2219] = v[1277] + v[1284] + v[1291] + v[1293] * v[1607] + v[1292] * v[1611] + v[1281] * v[1613] + v[1168] * v[1772]
		+ v[1167] * v[1777] + v[1166] * v[1782];
	v[1366] = (-(alphaB[0] * v[1353]) + alphaB[1] * v[1354] - alphaB[2] * v[1355] + 4e0*v[2219] * v[258] + v[1352] * v[513]
		+ 2e0*v[1349] * v[518] + 2e0*v[1347] * v[524] + 2e0*v[1340] * v[528] + v[1338] * v[532] + 2e0*v[1335] * v[537]
		+ 2e0*v[1334] * v[541] + 2e0*v[1332] * v[545] + v[1329] * v[550]) / 2e0;
	v[1414] = -v[1360] + v[1301] * v[1531] + v[1297] * v[1533] + v[1286] * v[1535] - 4e0*v[1366] * v[977];
	v[2285] = v[1414] - (v[1352] * v[258]) / 2e0;
	v[2283] = -v[1356] + v[1360] + v[2285];
	v[2278] = -v[1356] + v[1414];
	v[1225] = dGAp[2][0] * v[1156] + dGAp[2][1] * v[1157] + v[1225] + GAp[2] * v[1265];
	v[1226] = dGAp[1][0] * v[1156] + dGAp[1][1] * v[1157] + v[1226] + GAp[1] * v[1265];
	v[1227] = dGAp[0][0] * v[1156] + dGAp[0][1] * v[1157] + v[1227] + GAp[0] * v[1265];
	v[1367] = QAi[2][2] * v[1225] + QAi[2][1] * v[1226] + QAi[2][0] * v[1227] + v[2978] * v[318];
	v[1370] = QAi[1][2] * v[1225] + QAi[1][1] * v[1226] + QAi[1][0] * v[1227] + v[1320] * v[2712] + v[2806] * v[316];
	v[1372] = QAi[0][2] * v[1225] + QAi[0][1] * v[1226] + QAi[0][0] * v[1227] + v[2970] * v[310];
	v[1217] = dGAp[2][0] * v[1163] + dGAp[2][1] * v[1165] + v[1217] + GAp[2] * v[1267];
	v[1218] = dGAp[1][0] * v[1163] + dGAp[1][1] * v[1165] + v[1218] + GAp[1] * v[1267];
	v[1219] = dGAp[0][0] * v[1163] + dGAp[0][1] * v[1165] + v[1219] + GAp[0] * v[1267];
	v[1373] = QAi[2][2] * v[1217] + QAi[2][1] * v[1218] + QAi[2][0] * v[1219] + v[1319] * v[2710] + v[2807] * v[318];
	v[1376] = QAi[1][2] * v[1217] + QAi[1][1] * v[1218] + QAi[1][0] * v[1219] + v[2974] * v[316];
	v[1378] = QAi[0][2] * v[1217] + QAi[0][1] * v[1218] + QAi[0][0] * v[1219] + v[2969] * v[310];
	v[1209] = dGAp[2][0] * v[1174] + dGAp[2][1] * v[1175] + v[1209] + GAp[2] * v[1268];
	v[1210] = dGAp[1][0] * v[1174] + dGAp[1][1] * v[1175] + v[1210] + GAp[1] * v[1268];
	v[1211] = dGAp[0][0] * v[1174] + dGAp[0][1] * v[1175] + v[1211] + GAp[0] * v[1268];
	v[1385] = QAi[2][2] * v[1209] + QAi[2][1] * v[1210] + QAi[2][0] * v[1211] + v[1308] * v[1418] + v[2977] * v[318];
	v[1387] = QAi[1][2] * v[1209] + QAi[1][1] * v[1210] + QAi[1][0] * v[1211] + v[2808] * v[316];
	v[1390] = QAi[0][2] * v[1209] + QAi[0][1] * v[1210] + QAi[0][0] * v[1211] + v[2968] * v[310];
	v[1391] = -(v[1323] * v[418]) + 2e0*v[1318] * v[421] + v[1327] * v[424];
	v[2950] = -v[1391] / 2e0;
	v[1392] = 2e0*v[1311] * v[418] - v[1323] * v[421] - v[1312] * v[424];
	v[2949] = v[1392] / 2e0;
	v[1393] = -(v[1312] * v[418]) + v[1327] * v[421] + 2e0*v[1304] * v[424];
	v[2948] = -v[1393] / 2e0;
	v[4968] = 0e0;
	v[4969] = 0e0;
	v[4970] = 0e0;
	v[4971] = v[2950];
	v[4972] = v[2949];
	v[4973] = v[2948];
	v[4974] = 0e0;
	v[4975] = 0e0;
	v[4976] = 0e0;
	v[4977] = 0e0;
	v[4978] = 0e0;
	v[4979] = 0e0;
	v[1394] = (v[1367] * v[218]) / 2e0;
	v[1395] = v[1370] * v[218];
	v[1396] = v[1372] * v[218];
	v[1397] = v[1373] * v[218];
	v[1398] = (v[1376] * v[218]) / 2e0;
	v[1399] = v[1378] * v[218];
	v[1400] = v[1385] * v[218];
	v[1401] = v[1387] * v[218];
	v[1403] = 1e0 / Power(v[307], 2);
	v[2818] = -(v[1403] * v[316]);
	v[2817] = -(v[1403] * v[318]);
	v[2816] = -(v[1403] * v[310]);
	v[2815] = -(v[1403] * v[305]);
	v[2814] = -(v[1169] * v[1403]);
	v[2813] = -(v[1403] * v[306]);
	v[2812] = -(v[1403] * v[2709]);
	v[2263] = -(v[1403] * (v[2937] + v[223] * v[308]));
	v[2262] = -(v[1403] * (v[2933] + v[222] * v[308]));
	v[2261] = -(v[1403] * (v[2931] + v[2935]));
	v[2260] = -(v[1403] * (v[2936] + v[2938]));
	v[2259] = -(v[1403] * (v[2929] + v[2932]));
	v[2258] = -(v[1403] * (v[2930] + v[2934]));
	v[2257] = -(v[1403] * v[308]);
	v[2941] = v[2257] * v[310];
	v[2256] = -(v[1403] * v[311]);
	v[2940] = v[2256] * v[316];
	v[2939] = -(v[1403] * v[318] * v[320]);
	v[2254] = -(v[1403] * v[2809]);
	v[2253] = -(v[1403] * v[2810]);
	v[2252] = -(v[1403] * v[2811]);
	v[2251] = v[1320] * v[2815];
	v[2250] = v[1319] * v[2813];
	v[2249] = v[1308] * v[2812];
	v[2168] = v[1170] * v[2812];
	v[2167] = v[2711] * v[2814];
	v[2165] = v[1171] * v[2257];
	v[3165] = v[2249] + (v[2165] + v[2167])*v[318];
	v[2915] = v[2165] + v[2168];
	v[3166] = v[2167] + v[2915];
	v[2161] = v[1171] * v[2813];
	v[2158] = v[1170] * v[2256];
	v[3164] = v[2158] + v[2161];
	v[2157] = v[2713] * v[2814];
	v[2914] = v[2157] + v[2158];
	v[3163] = v[2161] + v[2914];
	v[3162] = v[2250] + v[2914] * v[318];
	v[2153] = v[1171] * v[2815];
	v[2151] = -(v[1170] * v[1403] * v[317]);
	v[2150] = v[2814] * v[320];
	v[3161] = v[2150] + v[2153];
	v[2913] = v[2150] + v[2151];
	v[3160] = v[2251] + v[2913] * v[316];
	v[3159] = v[2153] + v[2913];
	v[2071] = v[230] * v[2816];
	v[2068] = v[225] * v[2816];
	v[2065] = v[223] * v[2817];
	v[2064] = v[222] * v[2818];
	v[2059] = v[228] * v[2817];
	v[2058] = v[227] * v[2818];
	v[3149] = v[2058] + v[2068];
	v[2907] = v[2058] + v[2059];
	v[3147] = v[2068] + v[2907];
	v[2056] = v[221] * v[2816];
	v[3145] = v[2056] + v[2065];
	v[2908] = v[2056] + v[2064];
	v[3148] = v[2065] + v[2908];
	v[2052] = v[233] * v[2817];
	v[3150] = v[2052] + v[2071];
	v[2051] = v[232] * v[2818];
	v[2905] = v[2051] + v[2052];
	v[3146] = v[2071] + v[2905];
	v[1767] = v[2939] + v[2259] * v[310] + v[2258] * v[316];
	v[1762] = v[2940] + v[2261] * v[310] + v[2260] * v[318];
	v[1757] = v[2941] + v[2262] * v[316] + v[2263] * v[318];
	v[1585] = v[223] * v[2812];
	v[1583] = v[228] * v[2813];
	v[1579] = v[232] * v[2815];
	v[2248] = v[1304] + v[1311] + v[1318] + v[1320] * v[1579] + v[1319] * v[1583] + v[1308] * v[1585] + v[1171] * v[1757]
		+ v[1170] * v[1762] + v[1169] * v[1767];
	v[1404] = (-(alphaA[0] * v[1391]) + alphaA[1] * v[1392] - alphaA[2] * v[1393] + 4e0*v[218] * v[2248] + v[1390] * v[414]
		+ 2e0*v[1387] * v[419] + 2e0*v[1385] * v[425] + 2e0*v[1378] * v[429] + v[1376] * v[433] + 2e0*v[1373] * v[438]
		+ 2e0*v[1372] * v[442] + 2e0*v[1370] * v[446] + v[1367] * v[451]) / 2e0;
	v[1412] = -v[1398] + v[1328] * v[1480] + v[1324] * v[1482] + v[1313] * v[1484] - 4e0*v[1404] * v[962];
	v[2277] = v[1412] - (v[1390] * v[218]) / 2e0;
	v[2275] = -v[1394] + v[1398] + v[2277];
	v[2270] = -v[1394] + v[1412];
	v[2280] = (v[1358] + v[1362]) / 2e0;
	v[1406] = v[1357] + v[1359];
	v[2284] = v[1406] / 2e0;
	v[1407] = v[1361] + v[1363];
	v[2279] = v[1407] / 2e0;
	v[2272] = (v[1396] + v[1400]) / 2e0;
	v[1409] = v[1395] + v[1397];
	v[2276] = v[1409] / 2e0;
	v[1410] = v[1399] + v[1401];
	v[4208] = v[1097] + v[1208];
	v[4209] = v[1101] + v[1216];
	v[4210] = v[1105] + v[1224];
	v[4211] = v[1395] - v[1397] + 2e0*alphaA[0] * v[2270] + alphaA[2] * v[2272] + v[218] * v[2950] + v[1410] * v[400]
		- v[1106] * v[779] - v[1107] * v[780] - v[1108] * v[781] + v[2947] * v[962];
	v[4212] = -v[1396] + v[1400] + (alphaA[2] * v[1409]) / 2e0 + (alphaA[0] * v[1410]) / 2e0 + 2e0*alphaA[1] * v[2275]
		+ v[218] * v[2949] - v[1110] * v[779] - v[1111] * v[780] - v[1112] * v[781] + v[2946] * v[962];
	v[4213] = v[1399] - v[1401] + alphaA[0] * v[2272] + 2e0*alphaA[2] * v[2277] + v[218] * v[2948] + v[1409] * v[400]
		- v[1114] * v[779] - v[1115] * v[780] - v[1116] * v[781] + v[2945] * v[962];
	v[4214] = -v[1097] + v[1204];
	v[4215] = -v[1101] + v[1212];
	v[4216] = -v[1105] + v[1220];
	v[4217] = v[1357] - v[1359] + 2e0*alphaB[0] * v[2278] + alphaB[2] * v[2280] + v[258] * v[2956] + v[1407] * v[406]
		- v[1130] * v[779] - v[1131] * v[780] - v[1132] * v[781] + v[2953] * v[977];
	v[4218] = -v[1358] + v[1362] + (alphaB[2] * v[1406]) / 2e0 + (alphaB[0] * v[1407]) / 2e0 + 2e0*alphaB[1] * v[2283]
		+ v[258] * v[2955] - v[1134] * v[779] - v[1135] * v[780] - v[1136] * v[781] + v[2952] * v[977];
	v[4219] = v[1361] - v[1363] + alphaB[0] * v[2280] + 2e0*alphaB[2] * v[2285] + v[258] * v[2954] + v[1406] * v[406]
		- v[1138] * v[779] - v[1139] * v[780] - v[1140] * v[781] + v[2951] * v[977];
	v[2271] = v[1410] / 2e0;
	for (i1143 = 1; i1143 <= 12; i1143++) {
		i2850 = (i1143 == 9 ? 1 : 0);
		i2849 = (i1143 == 8 ? 1 : 0);
		i2848 = (i1143 == 7 ? 1 : 0);
		i2847 = (i1143 == 3 ? 1 : 0);
		i2846 = (i1143 == 2 ? 1 : 0);
		i2845 = (i1143 == 1 ? 1 : 0);
		i2835 = i2845 - i2848;
		i2834 = i2846 - i2849;
		i2833 = i2847 - i2850;
		i2824 = (i1143 == 11 ? 1 : 0);
		v[2894] = (*a4)*i2824;
		i2823 = (i1143 == 10 ? 1 : 0);
		v[2895] = (*a4)*i2823;
		i2822 = (i1143 == 12 ? 1 : 0);
		v[2896] = (*a4)*i2822;
		i2821 = (i1143 == 5 ? 1 : 0);
		v[2902] = (*a4)*i2821;
		i2820 = (i1143 == 4 ? 1 : 0);
		v[2903] = (*a4)*i2820;
		i2819 = (i1143 == 6 ? 1 : 0);
		v[2904] = (*a4)*i2819;
		v[1432] = v[3783 + i1143];
		v[1434] = v[3807 + i1143];
		v[1436] = v[3795 + i1143];
		v[1438] = v[3831 + i1143];
		v[1440] = v[3855 + i1143];
		v[1442] = v[3843 + i1143];
		v[1444] = v[3771 + i1143];
		v[2273] = -8e0*v[1411] * v[1444];
		v[1477] = 2e0*v[1444] * v[962];
		v[1497] = -4e0*v[1477] * v[218];
		v[2825] = v[1497] * v[307];
		v[1445] = v[4235 + i1143];
		v[1446] = v[4259 + i1143];
		v[1448] = v[3819 + i1143];
		v[2281] = -8e0*v[1413] * v[1448];
		v[1528] = 2e0*v[1448] * v[977];
		v[1548] = -4e0*v[1528] * v[258];
		v[2826] = v[1548] * v[333];
		v[1449] = v[4307 + i1143];
		v[1450] = v[4331 + i1143];
		v[1463] = -i2819 + v[1432];
		v[1464] = i2819 + v[1432];
		v[1465] = -i2820 + v[1434];
		v[1466] = i2820 + v[1434];
		v[1467] = i2821 + v[1436];
		v[1468] = -i2821 + v[1436];
		v[1469] = -i2822 + v[1438];
		v[1470] = i2822 + v[1438];
		v[1471] = -i2823 + v[1440];
		v[1472] = i2823 + v[1440];
		v[1473] = i2824 + v[1442];
		v[1474] = -i2824 + v[1442];
		v[1475] = 2e0*alphaA[1] * i2821 - v[1444];
		v[1476] = alphaA[2] * v[1477] - (i2819*v[218]) / 2e0;
		v[1478] = -(alphaA[1] * v[1477]) + (i2821*v[218]) / 2e0;
		v[1479] = alphaA[0] * v[1477] - (i2820*v[218]) / 2e0;
		v[1481] = v[1444] * v[1480] - 8e0*i2821*v[962];
		v[1483] = v[1444] * v[1482] + 8e0*i2819*v[962];
		v[1485] = v[1444] * v[1484] + 8e0*i2820*v[962];
		v[1486] = (v[1445] * v[218]) / 2e0 - v[1477] * v[414];
		v[1504] = v[1486] * v[310];
		v[1487] = v[1463] * v[218] - 2e0*v[1477] * v[419];
		v[1488] = v[1467] * v[218] - 2e0*v[1477] * v[425];
		v[1505] = v[1488] * v[318];
		v[1489] = v[1464] * v[218] - 2e0*v[1477] * v[429];
		v[1490] = (v[1475] * v[218]) / 2e0 - v[1477] * v[433];
		v[1512] = v[1490] * v[316];
		v[1491] = v[1465] * v[218] - 2e0*v[1477] * v[438];
		v[1513] = v[1491] * v[318];
		v[1492] = v[1468] * v[218] - 2e0*v[1477] * v[442];
		v[1521] = v[1492] * v[310];
		v[1493] = v[1466] * v[218] - 2e0*v[1477] * v[446];
		v[1494] = (v[1446] * v[218]) / 2e0 - v[1477] * v[451];
		v[1519] = v[1494] * v[318];
		v[1495] = v[1497] + 2e0*v[1476] * v[424];
		v[1496] = -(v[1476] * v[418]) - v[1478] * v[424];
		v[1498] = v[1497] + 2e0*v[1478] * v[418];
		v[1601] = v[1170] * v[1498];
		v[1499] = v[1476] * v[421] + v[1479] * v[424];
		v[1500] = -(v[1479] * v[418]) - v[1478] * v[421];
		v[1501] = v[1497] + 2e0*v[1479] * v[421];
		v[1597] = v[1171] * v[1501];
		v[1502] = v[1504] + v[1487] * v[316];
		v[1503] = v[1502] + v[1505] + v[2903];
		v[1506] = v[1504] + v[1505];
		v[1507] = QAi[0][0] * v[1486] + QAi[1][0] * v[1487] + QAi[2][0] * v[1488];
		v[1508] = QAi[0][1] * v[1486] + QAi[1][1] * v[1487] + QAi[2][1] * v[1488];
		v[1509] = QAi[0][2] * v[1486] + QAi[1][2] * v[1487] + QAi[2][2] * v[1488];
		v[1510] = v[1512] + v[1489] * v[310];
		v[1511] = v[1510] + v[1513] + v[2902];
		v[1514] = v[1512] + v[1513];
		v[1515] = QAi[0][0] * v[1489] + QAi[1][0] * v[1490] + QAi[2][0] * v[1491];
		v[1516] = QAi[0][1] * v[1489] + QAi[1][1] * v[1490] + QAi[2][1] * v[1491];
		v[1517] = QAi[0][2] * v[1489] + QAi[1][2] * v[1490] + QAi[2][2] * v[1491];
		v[1518] = v[1519] + v[1521];
		v[1520] = v[1519] + v[1493] * v[316];
		v[1522] = v[1520] + v[1521] + v[2904];
		v[1523] = QAi[0][0] * v[1492] + QAi[1][0] * v[1493] + QAi[2][0] * v[1494];
		v[2839] = v[1145] * v[1507] + v[1146] * v[1515] + v[1147] * v[1523];
		v[1524] = QAi[0][1] * v[1492] + QAi[1][1] * v[1493] + QAi[2][1] * v[1494];
		v[2840] = v[1145] * v[1508] + v[1146] * v[1516] + v[1147] * v[1524];
		v[1525] = QAi[0][2] * v[1492] + QAi[1][2] * v[1493] + QAi[2][2] * v[1494];
		v[2841] = v[1145] * v[1509] + v[1146] * v[1517] + v[1147] * v[1525];
		v[1526] = 2e0*alphaB[1] * i2824 - v[1448];
		v[1527] = alphaB[2] * v[1528] - (i2822*v[258]) / 2e0;
		v[1529] = -(alphaB[1] * v[1528]) + (i2824*v[258]) / 2e0;
		v[1530] = alphaB[0] * v[1528] - (i2823*v[258]) / 2e0;
		v[1532] = v[1448] * v[1531] - 8e0*i2824*v[977];
		v[1534] = v[1448] * v[1533] + 8e0*i2822*v[977];
		v[1536] = v[1448] * v[1535] + 8e0*i2823*v[977];
		v[1537] = (v[1449] * v[258]) / 2e0 - v[1528] * v[513];
		v[1555] = v[1537] * v[336];
		v[1538] = v[1469] * v[258] - 2e0*v[1528] * v[518];
		v[1539] = v[1473] * v[258] - 2e0*v[1528] * v[524];
		v[1556] = v[1539] * v[344];
		v[1540] = v[1470] * v[258] - 2e0*v[1528] * v[528];
		v[1541] = (v[1526] * v[258]) / 2e0 - v[1528] * v[532];
		v[1563] = v[1541] * v[342];
		v[1542] = v[1471] * v[258] - 2e0*v[1528] * v[537];
		v[1564] = v[1542] * v[344];
		v[1543] = v[1474] * v[258] - 2e0*v[1528] * v[541];
		v[1572] = v[1543] * v[336];
		v[1544] = v[1472] * v[258] - 2e0*v[1528] * v[545];
		v[1545] = (v[1450] * v[258]) / 2e0 - v[1528] * v[550];
		v[1570] = v[1545] * v[344];
		v[1546] = v[1548] + 2e0*v[1527] * v[523];
		v[1547] = -(v[1527] * v[517]) - v[1529] * v[523];
		v[1549] = v[1548] + 2e0*v[1529] * v[517];
		v[1629] = v[1167] * v[1549];
		v[1550] = v[1527] * v[520] + v[1530] * v[523];
		v[1551] = -(v[1530] * v[517]) - v[1529] * v[520];
		v[1552] = v[1548] + 2e0*v[1530] * v[520];
		v[1625] = v[1168] * v[1552];
		v[1553] = v[1555] + v[1538] * v[342];
		v[1554] = v[1553] + v[1556] + v[2895];
		v[1557] = v[1555] + v[1556];
		v[1558] = QBi[0][0] * v[1537] + QBi[1][0] * v[1538] + QBi[2][0] * v[1539];
		v[1559] = QBi[0][1] * v[1537] + QBi[1][1] * v[1538] + QBi[2][1] * v[1539];
		v[1560] = QBi[0][2] * v[1537] + QBi[1][2] * v[1538] + QBi[2][2] * v[1539];
		v[1561] = v[1563] + v[1540] * v[336];
		v[1562] = v[1561] + v[1564] + v[2894];
		v[1565] = v[1563] + v[1564];
		v[1566] = QBi[0][0] * v[1540] + QBi[1][0] * v[1541] + QBi[2][0] * v[1542];
		v[1567] = QBi[0][1] * v[1540] + QBi[1][1] * v[1541] + QBi[2][1] * v[1542];
		v[1568] = QBi[0][2] * v[1540] + QBi[1][2] * v[1541] + QBi[2][2] * v[1542];
		v[1569] = v[1570] + v[1572];
		v[1571] = v[1570] + v[1544] * v[342];
		v[1573] = v[1571] + v[1572] + v[2896];
		v[1574] = QBi[0][0] * v[1543] + QBi[1][0] * v[1544] + QBi[2][0] * v[1545];
		v[2842] = v[1145] * v[1558] + v[1146] * v[1566] + v[1147] * v[1574];
		v[1575] = QBi[0][1] * v[1543] + QBi[1][1] * v[1544] + QBi[2][1] * v[1545];
		v[2843] = v[1145] * v[1559] + v[1146] * v[1567] + v[1147] * v[1575];
		v[1576] = QBi[0][2] * v[1543] + QBi[1][2] * v[1544] + QBi[2][2] * v[1545];
		v[2844] = v[1145] * v[1560] + v[1146] * v[1568] + v[1147] * v[1576];
		v[1577] = v[1481] + v[1499];
		v[1595] = v[1171] * v[1577];
		v[1578] = -v[1481] + v[1499];
		v[1580] = (v[1577] * v[232] + v[1579] * v[2825] + v[1493] * v[305]) / v[307];
		v[1581] = v[1483] + v[1500];
		v[1603] = v[1171] * v[1581];
		v[1582] = -v[1483] + v[1500];
		v[1599] = v[1170] * v[1582];
		v[1584] = (v[1581] * v[228] + v[1583] * v[2825] + v[1491] * v[306]) / v[307];
		v[1586] = (v[1582] * v[223] + v[1488] * v[2709] + v[1585] * v[2825]) / v[307];
		v[1587] = v[1597] + v[1599];
		v[2975] = v[1587] / v[307];
		v[1588] = v[1485] + v[1496];
		v[1593] = v[1170] * v[1588];
		v[1589] = -v[1485] + v[1496];
		v[1590] = v[1601] + v[1603];
		v[2971] = v[1590] / v[307];
		v[1591] = v[1169] * v[1495] + v[1593] + v[1595];
		v[2979] = v[1591] / v[307];
		v[1594] = v[1591] - v[1595];
		v[1596] = v[1591] - v[1593];
		v[2972] = v[1596] / v[307];
		v[1598] = v[1169] * v[1578] + v[1597];
		v[1600] = v[1598] + v[1599];
		v[2973] = v[1600] / v[307];
		v[1602] = v[1169] * v[1589] + v[1601];
		v[1604] = v[1602] + v[1603];
		v[2976] = v[1604] / v[307];
		v[1605] = v[1532] + v[1550];
		v[1623] = v[1168] * v[1605];
		v[1606] = -v[1532] + v[1550];
		v[1608] = (v[1605] * v[272] + v[1607] * v[2826] + v[1544] * v[331]) / v[333];
		v[1609] = v[1534] + v[1551];
		v[1631] = v[1168] * v[1609];
		v[1610] = -v[1534] + v[1551];
		v[1627] = v[1167] * v[1610];
		v[1612] = (v[1609] * v[268] + v[1611] * v[2826] + v[1542] * v[332]) / v[333];
		v[1614] = (v[1610] * v[263] + v[1539] * v[2714] + v[1613] * v[2826]) / v[333];
		v[1615] = v[1625] + v[1627];
		v[2987] = v[1615] / v[333];
		v[1616] = v[1536] + v[1547];
		v[1621] = v[1167] * v[1616];
		v[1617] = -v[1536] + v[1547];
		v[1618] = v[1629] + v[1631];
		v[2983] = v[1618] / v[333];
		v[1619] = v[1166] * v[1546] + v[1621] + v[1623];
		v[2991] = v[1619] / v[333];
		v[1622] = v[1619] - v[1623];
		v[1624] = v[1619] - v[1621];
		v[2984] = v[1624] / v[333];
		v[1626] = v[1166] * v[1606] + v[1625];
		v[1628] = v[1626] + v[1627];
		v[2985] = v[1628] / v[333];
		v[1630] = v[1166] * v[1617] + v[1629];
		v[1632] = v[1630] + v[1631];
		v[2988] = v[1632] / v[333];
		v[1634] = i2835 + GAp[0] * v[1507] + GAp[1] * v[1508] + GAp[2] * v[1509] - GBp[0] * v[1558] - GBp[1] * v[1559]
			- GBp[2] * v[1560];
		v[1637] = i2834 + GAp[0] * v[1515] + GAp[1] * v[1516] + GAp[2] * v[1517] - GBp[0] * v[1566] - GBp[1] * v[1567]
			- GBp[2] * v[1568];
		v[1640] = i2833 + GAp[0] * v[1523] + GAp[1] * v[1524] + GAp[2] * v[1525] - GBp[0] * v[1574] - GBp[1] * v[1575]
			- GBp[2] * v[1576];
		v[2892] = v[1634] * v[356] + v[1637] * v[357] + v[1640] * v[358];
		v[1642] = v[2892] / v[1264];
		v[2827] = v[1642] * v[369];
		v[1652] = v[2827] * v[358] + v[1640] * v[370];
		v[1648] = v[2827] * v[357] + v[1637] * v[370];
		v[1644] = v[2827] * v[356] + v[1634] * v[370];
		v[1645] = 2e0*v[1201] * v[1644];
		v[1646] = v[1644];
		v[1649] = 2e0*v[1193] * v[1648];
		v[1650] = v[1648];
		v[1653] = 2e0*v[1192] * v[1652];
		v[1654] = v[1652];
		v[1655] = 0e0;
		v[1656] = 0e0;
		v[1657] = 0e0;
		v[1658] = 0e0;
		v[1659] = 0e0;
		v[1660] = 0e0;
		v[1661] = 0e0;
		v[1662] = 0e0;
		v[1663] = 0e0;
		v[1664] = 0e0;
		v[1665] = 0e0;
		v[1666] = 0e0;
		v[1667] = 0e0;
		v[1668] = 0e0;
		v[1669] = 0e0;
		v[1670] = 0e0;
		v[1671] = 0e0;
		v[1672] = 0e0;
		v[1673] = 0e0;
		v[1674] = 0e0;
		v[1675] = 0e0;
		v[1676] = 0e0;
		v[1677] = 0e0;
		v[1678] = 0e0;
		v[1679] = 0e0;
		v[1680] = 0e0;
		v[1681] = 0e0;
		v[1682] = 0e0;
		v[1683] = 0e0;
		v[1684] = 0e0;
		v[1685] = 0e0;
		v[1686] = 0e0;
		b1687 = b709;
		if (b1687) {
			v[1690] = v[1654];
			v[1688] = v[1650];
			v[1689] = v[1654] * v[381] - v[1688] * v[382];
			v[1691] = -(v[1690] * v[380]) + v[1646] * v[382];
			v[1692] = v[1688] * v[380] - v[1646] * v[381];
			v[2831] = v[1255] * v[1689] + v[1250] * v[1691] + v[1245] * v[1692];
			v[2829] = v[1689] * v[711] + v[1691] * v[712] + v[1692] * v[713];
			v[1693] = v[2829] / v[714];
			v[2830] = v[1693] * v[3063];
			v[1703] = v[1693] * v[2828];
			v[1658] = -(v[1258] * v[1694] * v[2829]);
			v[1695] = v[1689] * v[2762] + v[2830] * v[711];
			v[1714] = 2e0*v[1695] * v[723];
			v[1697] = v[1691] * v[2762] + v[2830] * v[712];
			v[1711] = 2e0*v[1697] * v[724];
			v[1698] = v[1692] * v[2762] + v[2830] * v[713];
			v[1712] = 2e0*v[1698] * v[725];
			v[1661] = v[1703] * v[2790] + v[2831] * v[716];
			v[1660] = v[1693] * v[1701] * v[720];
			v[1659] = v[1693] * v[2790] * v[721] + v[2831] * v[722];
			v[1685] = 2e0*v[1703] * v[1704] * v[3064];
			v[1686] = v[1257] * v[1693] * v[1699] * v[1704] * v[3065];
			v[1657] = v[1245] * v[2830] + v[1692] * v[2832];
			v[1656] = v[1250] * v[2830] + v[1691] * v[2832];
			v[1655] = v[1255] * v[2830] + v[1689] * v[2832];
			v[1706] = (v[1697] * v[723] + v[1695] * v[724]) / 2e0;
			v[1707] = v[1711] + v[1714];
			v[1708] = v[1707] + v[1712];
			v[1709] = (v[1698] * v[723] + v[1695] * v[725]) / 2e0;
			v[1710] = (v[1698] * v[724] + v[1697] * v[725]) / 2e0;
			v[1713] = v[1711] + v[1712];
			v[1715] = v[1712] + v[1714];
			v[1664] = (v[1230] * v[1695] + v[1229] * v[1697] + 4e0*v[1698] * v[1716]) / 2e0;
			v[1663] = (v[1231] * v[1695] + v[1229] * v[1698] + 4e0*v[1697] * v[1717]) / 2e0;
			v[1662] = (v[1231] * v[1697] + v[1230] * v[1698] + 4e0*v[1695] * v[1718]) / 2e0;
			v[1719] = -4e0*v[1708] * v[1962];
			v[1684] = 8e0*v[1241] * v[1708] * v[3066];
			v[1675] = (v[1182] * v[1719]) / 2e0;
			v[1679] = (v[1186] * v[1719]) / 2e0;
			v[1677] = v[1184] * v[1719];
			v[1681] = v[1188] * v[1719];
			v[1680] = v[1187] * v[1719];
			v[1682] = v[1189] * v[1719];
			v[1676] = v[1183] * v[1719];
			v[1678] = v[1185] * v[1719];
			v[1683] = (v[1190] * v[1719]) / 2e0;
			v[1720] = -v[1698] + v[1706];
			v[1721] = v[1698] + v[1706];
			v[1722] = v[1697] + v[1709];
			v[1723] = -v[1697] + v[1709];
			v[1724] = -v[1695] + v[1710];
			v[1725] = v[1695] + v[1710];
			v[1667] = v[1233] * v[1719] + v[1720] * v[726];
			v[1669] = v[1235] * v[1719] + v[1721] * v[726];
			v[1668] = v[1234] * v[1719] + v[1722] * v[726];
			v[1672] = v[1238] * v[1719] + v[1723] * v[726];
			v[1666] = (v[1232] * v[1719] - v[1713] * v[726]) / 2e0;
			v[1671] = v[1237] * v[1719] + v[1724] * v[726];
			v[1673] = v[1239] * v[1719] + v[1725] * v[726];
			v[1670] = (v[1236] * v[1719] - v[1715] * v[726]) / 2e0;
			v[1674] = (v[1240] * v[1719] - v[1707] * v[726]) / 2e0;
			v[1665] = -(v[1190] * v[1707]) / 2e0 - (v[1182] * v[1713]) / 2e0 - (v[1186] * v[1715]) / 2e0 + v[1183] * v[1720]
				+ v[1185] * v[1721] + v[1184] * v[1722] + v[1188] * v[1723] + v[1187] * v[1724] + v[1189] * v[1725];
		}
		else {
		};
		v[1726] = 0e0;
		v[1727] = 0e0;
		v[1728] = 0e0;
		v[1729] = 0e0;
		v[1730] = 0e0;
		v[1731] = 0e0;
		b1732 = (*previouscontact);
		if (b1732) {
			v[1750] = v[1654];
			v[1744] = v[1650];
			v[1742] = v[1646];
			v[1740] = v[1667];
			v[1739] = v[1668];
			v[1737] = v[1670];
			v[1736] = v[1671];
			v[1734] = v[1673];
			v[1733] = v[1674];
			v[1748] = v[1178] * v[1742];
			v[1747] = v[1176] * v[1654];
			v[1745] = v[1177] * v[1650];
			v[1674] = 0e0;
			v[1673] = 0e0;
			v[1735] = i2833 + GAi[0] * v[1523] + GAi[1] * v[1524] + GAi[2] * v[1525] - GBi[0] * v[1574] - GBi[1] * v[1575]
				- GBi[2] * v[1576] + gti[0] * v[1672] + gti[2] * v[1733] + gti[1] * v[1734];
			v[2836] = -(v[1735] * v[373]) - v[1750] * v[746];
			v[1672] = 0e0;
			v[1671] = 0e0;
			v[1670] = 0e0;
			v[1738] = i2834 + GAi[0] * v[1515] + GAi[1] * v[1516] + GAi[2] * v[1517] - GBi[0] * v[1566] - GBi[1] * v[1567]
				- GBi[2] * v[1568] + gti[0] * v[1669] + gti[2] * v[1736] + gti[1] * v[1737];
			v[2838] = -(v[1738] * v[372]) - v[1744] * v[745];
			v[1669] = 0e0;
			v[1668] = 0e0;
			v[1667] = 0e0;
			v[1741] = i2835 + GAi[0] * v[1507] + GAi[1] * v[1508] + GAi[2] * v[1509] - GBi[0] * v[1558] - GBi[1] * v[1559]
				- GBi[2] * v[1560] + gti[0] * v[1666] + gti[2] * v[1739] + gti[1] * v[1740];
			v[2837] = -(v[1741] * v[371]) - v[1742] * v[744];
			v[1666] = 0e0;
			v[1646] = 0e0;
			v[1743] = v[1745] + v[1748];
			v[1650] = 0e0;
			v[1746] = v[1745] + v[1747];
			v[1749] = v[1747] + v[1748];
			v[1654] = 0e0;
			v[1731] = v[1176] * v[1735];
			v[1730] = v[1177] * v[1738];
			v[1729] = v[1178] * v[1741];
			v[1726] = v[1742] * v[1751] - (2e0*v[1178] * v[1644] + v[1746])*v[371];
			v[1645] = v[1645] + v[1741] * v[1751] + v[1178] * (v[2836] + v[2838]) - v[1746] * v[744];
			v[1727] = v[1744] * v[1753] - (2e0*v[1177] * v[1648] + v[1749])*v[372];
			v[1649] = v[1649] + v[1738] * v[1753] + v[1177] * (v[2836] + v[2837]) - v[1749] * v[745];
			v[1728] = v[1750] * v[1756] - (2e0*v[1176] * v[1652] + v[1743])*v[373];
			v[1653] = v[1653] + v[1735] * v[1756] + v[1176] * (v[2837] + v[2838]) - v[1743] * v[746];
		}
		else {
		};
		v[1761] = v[1497] * v[1757] + v[1501] * v[1758] + v[1577] * v[1759] + v[1581] * v[1760] + v[1503] * v[2293]
			+ v[1510] * v[2710] + v[1518] * v[2712] + v[1580] * v[316] + v[1584] * v[318] + (*a4)*v[4523 + i1143];
		v[1766] = v[1418] * v[1502] + v[1417] * v[1520] + v[1497] * v[1762] + v[1498] * v[1763] + v[1582] * v[1764]
			+ v[1588] * v[1765] + v[1511] * v[2294] + v[1580] * v[310] + v[1586] * v[318] + (*a4)*v[4535 + i1143];
		v[1771] = v[1420] * v[1506] + v[1419] * v[1514] + v[1497] * v[1767] + v[1495] * v[1768] + v[1578] * v[1769]
			+ v[1589] * v[1770] + v[1522] * v[2295] + v[1584] * v[310] + v[1586] * v[316] + (*a4)*v[4547 + i1143];
		v[1776] = v[1548] * v[1772] + v[1552] * v[1773] + v[1605] * v[1774] + v[1609] * v[1775] + v[1554] * v[2296]
			+ v[1561] * v[2715] + v[1569] * v[2717] + v[1608] * v[342] + v[1612] * v[344] + (*a4)*v[4559 + i1143];
		v[1781] = v[1424] * v[1553] + v[1423] * v[1571] + v[1548] * v[1777] + v[1549] * v[1778] + v[1610] * v[1779]
			+ v[1616] * v[1780] + v[1562] * v[2297] + v[1608] * v[336] + v[1614] * v[344] + (*a4)*v[4571 + i1143];
		v[1786] = v[1426] * v[1557] + v[1425] * v[1565] + v[1548] * v[1782] + v[1546] * v[1783] + v[1606] * v[1784]
			+ v[1617] * v[1785] + v[1573] * v[2298] + v[1612] * v[336] + v[1614] * v[342] + (*a4)*v[4583 + i1143];
		v[1787] = (*ct)*(dGAp[0][1] * v[2839] + dGAp[1][1] * v[2840] + dGAp[2][1] * v[2841]);
		v[1788] = (*ct)*(dGAp[0][0] * v[2839] + dGAp[1][0] * v[2840] + dGAp[2][0] * v[2841]);
		v[1789] = (*ct)*(dGBp[0][1] * v[2842] + dGBp[1][1] * v[2843] + dGBp[2][1] * v[2844]);
		v[1790] = (*ct)*(dGBp[0][0] * v[2842] + dGBp[1][0] * v[2843] + dGBp[2][0] * v[2844]);
		v[1791] = -(i2845*v[1094]) + i2848 * v[1094] - i2846 * v[1098] + i2849 * v[1098] - i2847 * v[1102] + i2850 * v[1102]
			- i2820 * v[1106] - i2821 * v[1110] - i2819 * v[1114] - i2823 * v[1130] - i2824 * v[1134] - i2822 * v[1138];
		v[1806] = v[1791];
		v[1792] = -(i2845*v[1095]) + i2848 * v[1095] - i2846 * v[1099] + i2849 * v[1099] - i2847 * v[1103] + i2850 * v[1103]
			- i2820 * v[1107] - i2821 * v[1111] - i2819 * v[1115] - i2823 * v[1131] - i2824 * v[1135] - i2822 * v[1139];
		v[1804] = v[1792];
		v[1793] = -(i2845*v[1096]) + i2848 * v[1096] - i2846 * v[1100] + i2849 * v[1100] - i2847 * v[1104] + i2850 * v[1104]
			- i2820 * v[1108] - i2821 * v[1112] - i2819 * v[1116] - i2823 * v[1132] - i2824 * v[1136] - i2822 * v[1140];
		v[1802] = v[1793];
		v[1794] = 0e0;
		v[1795] = 0e0;
		v[1796] = 0e0;
		v[1797] = 0e0;
		v[1798] = 0e0;
		v[1799] = 0e0;
		b1800 = (*stick);
		if (b1800) {
			b1801 = b777;
			if (b1801) {
				v[1799] = v[1793];
				v[1793] = 0e0;
				v[1798] = v[1792];
				v[1792] = 0e0;
				v[1797] = v[1791];
				v[1791] = 0e0;
			}
			else {
				v[2851] = (*mud)*(v[1806] * v[773] + v[1804] * v[774] + v[1802] * v[775]);
				v[1793] = 0e0;
				v[1792] = 0e0;
				v[2853] = v[1817] * v[2851] * v[788];
				v[1791] = 0e0;
				v[2852] = v[1815] * v[2851] * v[787] * v[792];
				v[1799] = v[1802] * v[1814] + v[2852] * v[775];
				v[1798] = v[1804] * v[1814] + v[2852] * v[774];
				v[1797] = v[1806] * v[1814] + v[2852] * v[773];
				v[1796] = v[2853] * v[766];
				v[1795] = v[2853] * v[765];
				v[1794] = v[2853] * v[764];
			};
		}
		else {
			b1818 = b811;
			if (b1818) {
				v[1799] = v[1802];
				v[1793] = 0e0;
				v[1798] = v[1804];
				v[1792] = 0e0;
				v[1797] = v[1806];
				v[1791] = 0e0;
			}
			else {
				v[1824] = v[1806] * v[2854];
				v[1822] = v[1804] * v[2854];
				v[1821] = v[1802] * v[2854];
				v[2856] = (*mud)*v[1827] * (v[1806] * v[773] + v[1804] * v[774] + v[1802] * v[775])*v[818];
				v[2855] = v[1823] * (v[1824] * v[773] + v[1822] * v[774] + v[1821] * v[775])*v[817];
				v[1799] = v[2855] * v[775] + v[1821] * v[818];
				v[1798] = v[2855] * v[774] + v[1822] * v[818];
				v[1797] = v[2855] * v[773] + v[1824] * v[818];
				v[1796] = v[2856] * v[766];
				v[1795] = v[2856] * v[765];
				v[1794] = v[2856] * v[764];
			};
		};
		v[2869] = v[1794] * v[373];
		v[2867] = v[1794] * v[372];
		v[2865] = (*cn)*v[1794];
		v[1900] = -(v[2865] * v[749]);
		v[2864] = (*cn)*v[1795];
		v[2292] = -(v[1795] * v[1909]);
		v[2866] = v[1795] * v[372] + v[1796] * v[373];
		v[2863] = (*cn)*v[1796];
		v[1912] = -(v[2863] * v[753]);
		v[2870] = v[1912] - v[2292];
		v[2868] = v[1796] * v[1909] - v[2864] * v[751];
		v[2861] = -((*ct)*v[1797]);
		v[2859] = -((*ct)*v[1798]);
		v[2967] = (*ct)*(-(v[1102] * v[1797]) - v[1103] * v[1798] - v[1104] * v[1799]) - v[1912] + v[1794] * v[2290] + v[2292]
			- v[1788] * v[659] - v[1787] * v[674] + v[1790] * v[689] + v[1789] * v[704];
		v[2966] = (*ct)*(v[1098] * v[1797] + v[1099] * v[1798] + v[1100] * v[1799]) - v[1794] * v[2288] + v[2868]
			+ v[1788] * v[657] + v[1787] * v[672] - v[1790] * v[687] - v[1789] * v[702];
		v[2965] = (*ct)*(-(v[1094] * v[1797]) - v[1095] * v[1798] - v[1096] * v[1799]) - v[1900] + v[1795] * v[2288]
			+ v[1796] * v[2290] - v[1788] * v[655] - v[1787] * v[670] + v[1790] * v[685] + v[1789] * v[700];
		v[2857] = -((*ct)*v[1799]);
		v[1829] = -((*ct)*(v[1147] * v[1786] + v[1799] * v[355])) - i2822 * v[781];
		v[1830] = -((*ct)*(v[1147] * v[1781] + v[1799] * v[345])) - i2824 * v[781];
		v[1831] = -((*ct)*(v[1147] * v[1776] + v[1799] * v[335])) - i2823 * v[781];
		v[1836] = -((*ct)*(v[1147] * v[1771] + v[1799] * v[329])) - i2819 * v[781];
		v[1837] = -((*ct)*(v[1147] * v[1766] + v[1799] * v[319])) - i2821 * v[781];
		v[1838] = -((*ct)*(v[1147] * v[1761] + v[1799] * v[309])) - i2820 * v[781];
		v[2889] = -(i2847*v[2858]) + i2850 * v[2858] - v[2857] * v[300] + v[2857] * v[303];
		v[2890] = -(i2846*v[2858]) + i2849 * v[2858] - v[2857] * v[299] + v[2857] * v[302];
		v[2891] = -(i2845*v[2858]) + i2848 * v[2858] - v[2857] * v[298] + v[2857] * v[301];
		v[1842] = -((*ct)*(v[1146] * v[1786] + v[1798] * v[355])) - i2822 * v[780];
		v[1843] = -((*ct)*(v[1146] * v[1781] + v[1798] * v[345])) - i2824 * v[780];
		v[1844] = -((*ct)*(v[1146] * v[1776] + v[1798] * v[335])) - i2823 * v[780];
		v[1849] = -((*ct)*(v[1146] * v[1771] + v[1798] * v[329])) - i2819 * v[780];
		v[1850] = -((*ct)*(v[1146] * v[1766] + v[1798] * v[319])) - i2821 * v[780];
		v[1851] = -((*ct)*(v[1146] * v[1761] + v[1798] * v[309])) - i2820 * v[780];
		v[2886] = -(i2847*v[2860]) + i2850 * v[2860] - v[2859] * v[300] + v[2859] * v[303];
		v[2887] = -(i2846*v[2860]) + i2849 * v[2860] - v[2859] * v[299] + v[2859] * v[302];
		v[2888] = -(i2845*v[2860]) + i2848 * v[2860] - v[2859] * v[298] + v[2859] * v[301];
		v[1855] = -((*ct)*(v[1145] * v[1786] + v[1797] * v[355])) - i2822 * v[779];
		v[1856] = -((*ct)*(v[1145] * v[1781] + v[1797] * v[345])) - i2824 * v[779];
		v[1857] = -((*ct)*(v[1145] * v[1776] + v[1797] * v[335])) - i2823 * v[779];
		v[1862] = -((*ct)*(v[1145] * v[1771] + v[1797] * v[329])) - i2819 * v[779];
		v[1863] = -((*ct)*(v[1145] * v[1766] + v[1797] * v[319])) - i2821 * v[779];
		v[1864] = -((*ct)*(v[1145] * v[1761] + v[1797] * v[309])) - i2820 * v[779];
		v[2883] = -(i2847*v[2862]) + i2850 * v[2862] - v[2861] * v[300] + v[2861] * v[303];
		v[2884] = -(i2846*v[2862]) + i2849 * v[2862] - v[2861] * v[299] + v[2861] * v[302];
		v[2885] = -(i2845*v[2862]) + i2848 * v[2862] - v[2861] * v[298] + v[2861] * v[301];
		v[1871] = v[2474] * v[2863];
		v[1872] = v[2477] * v[2864];
		v[1653] = v[1653] + (*cn)*(v[1794] * v[2471] + v[1795] * v[2472] + v[1796] * v[2473]);
		v[1649] = v[1649] + (*cn)*(v[1794] * v[2476] + v[1795] * v[2478] + v[1796] * v[2479]);
		v[1881] = v[2481] * v[2865];
		v[1645] = v[1645] + (*cn)*(v[1794] * v[2482] + v[1795] * v[2483] + v[1796] * v[2484]);
		v[1899] = v[1898] * v[2866] + v[1900] * v[345];
		v[1902] = v[1901] * v[2866] + v[1900] * v[355];
		v[1904] = v[1903] * v[2866] + v[1900] * v[335];
		v[1905] = v[1903] * v[2867] + v[2868] * v[335];
		v[1908] = v[1901] * v[2867] + v[2868] * v[355];
		v[1910] = v[1901] * v[2869] + v[2870] * v[355];
		v[1911] = v[1898] * v[2867] + v[2868] * v[345];
		v[1913] = v[1903] * v[2869] + v[2870] * v[335];
		v[1914] = v[1898] * v[2869] + v[2870] * v[345];
		v[1916] = v[1915] * v[2866] - v[1900] * v[319];
		v[1918] = v[1917] * v[2866] - v[1900] * v[329];
		v[1920] = v[1919] * v[2866] - v[1900] * v[309];
		v[1921] = v[1919] * v[2867] - v[2868] * v[309];
		v[1922] = v[1917] * v[2867] - v[2868] * v[329];
		v[1923] = v[1917] * v[2869] - v[2870] * v[329];
		v[1924] = v[1915] * v[2867] - v[2868] * v[319];
		v[1925] = v[1919] * v[2869] - v[2870] * v[309];
		v[1926] = v[1915] * v[2869] - v[2870] * v[319];
		v[1927] = 0e0;
		v[1928] = 0e0;
		v[1929] = 0e0;
		v[1930] = 0e0;
		v[1931] = 0e0;
		v[1932] = 0e0;
		v[1933] = 0e0;
		v[1934] = 0e0;
		v[1935] = 0e0;
		b1936 = (*previouscontact);
		if (b1936) {
			v[1870] = (*epst)*v[1797];
			v[2873] = v[1870] * v[371];
			v[1869] = (*epst)*v[1798];
			v[2872] = v[1869] * v[372];
			v[2874] = v[2872] + v[2873];
			v[1868] = (*epst)*v[1799];
			v[2871] = v[1868] * v[373];
			v[2876] = v[2871] + v[2872];
			v[2875] = v[2871] + v[2873];
			v[1731] = v[1731] + v[1868] * v[746];
			v[1730] = v[1730] + v[1869] * v[745];
			v[1726] = v[1726] + v[1194] * v[1870] - v[2876] * v[371];
			v[1727] = v[1727] + v[1196] * v[1869] - v[2875] * v[372];
			v[1728] = v[1728] + v[1198] * v[1868] - v[2874] * v[373];
			v[1729] = v[1729] + v[1870] * v[744];
			v[1653] = v[1653] + v[1868] * v[2787] - v[2874] * v[746];
			v[1649] = v[1649] + v[1869] * v[2788] - v[2875] * v[745];
			v[1645] = v[1645] - v[1870] * v[2789] - v[2876] * v[744];
			v[1927] = gti[0] * v[1726];
			v[1928] = gti[1] * v[1726];
			v[1929] = gti[2] * v[1726];
			v[1937] = -v[1726];
			v[1938] = -(GBi[2] * v[1726]);
			v[1939] = -(GBi[1] * v[1726]);
			v[1940] = -(GBi[0] * v[1726]);
			v[1941] = v[1726];
			v[1942] = GAi[2] * v[1726];
			v[1943] = GAi[1] * v[1726];
			v[1944] = GAi[0] * v[1726];
			v[1930] = gti[0] * v[1727];
			v[1931] = gti[1] * v[1727];
			v[1932] = gti[2] * v[1727];
			v[1945] = -v[1727];
			v[1946] = -(GBi[2] * v[1727]);
			v[1947] = -(GBi[1] * v[1727]);
			v[1948] = -(GBi[0] * v[1727]);
			v[1949] = v[1727];
			v[1950] = GAi[2] * v[1727];
			v[1951] = GAi[1] * v[1727];
			v[1952] = GAi[0] * v[1727];
			v[1933] = gti[0] * v[1728];
			v[1934] = gti[1] * v[1728];
			v[1935] = gti[2] * v[1728];
			v[1953] = -v[1728];
			v[1954] = -(GBi[2] * v[1728]);
			v[1955] = -(GBi[1] * v[1728]);
			v[1956] = -(GBi[0] * v[1728]);
			v[1957] = v[1728];
			v[1958] = GAi[2] * v[1728];
			v[1959] = GAi[1] * v[1728];
			v[1960] = GAi[0] * v[1728];
			v[1881] = -v[1729] + v[1881];
			v[1872] = -v[1730] + v[1872];
			v[1871] = -v[1731] + v[1871];
		}
		else {
			v[1944] = 0e0;
			v[1943] = 0e0;
			v[1942] = 0e0;
			v[1952] = 0e0;
			v[1951] = 0e0;
			v[1950] = 0e0;
			v[1960] = 0e0;
			v[1959] = 0e0;
			v[1958] = 0e0;
			v[1941] = 0e0;
			v[1949] = 0e0;
			v[1957] = 0e0;
			v[1940] = 0e0;
			v[1939] = 0e0;
			v[1938] = 0e0;
			v[1948] = 0e0;
			v[1947] = 0e0;
			v[1946] = 0e0;
			v[1956] = 0e0;
			v[1955] = 0e0;
			v[1954] = 0e0;
			v[1937] = 0e0;
			v[1945] = 0e0;
			v[1953] = 0e0;
		};
		b1961 = b709;
		if (b1961) {
			v[1683] = v[1683] + (v[1935] * v[726]) / 2e0;
			v[1682] = v[1682] + v[1934] * v[726];
			v[1681] = v[1681] + v[1933] * v[726];
			v[1680] = v[1680] + v[1932] * v[726];
			v[1679] = v[1679] + (v[1931] * v[726]) / 2e0;
			v[1678] = v[1678] + v[1930] * v[726];
			v[1677] = v[1677] + v[1929] * v[726];
			v[1676] = v[1676] + v[1928] * v[726];
			v[1665] = v[1665] + (v[1232] * v[1927]) / 2e0 + v[1233] * v[1928] + v[1234] * v[1929] + v[1235] * v[1930] +
				(v[1236] * v[1931]) / 2e0 + v[1237] * v[1932] + v[1238] * v[1933] + v[1239] * v[1934] + (v[1240] * v[1935]) / 2e0;
			v[1675] = v[1675] + (v[1927] * v[726]) / 2e0;
			v[1684] = v[1684] - 4e0*v[1665] * v[1962];
			v[2877] = v[1675] - v[1684];
			v[2879] = (v[1677] + v[1681]) / 2e0;
			v[2878] = (v[1680] + v[1682]) / 2e0;
			v[1664] = v[1664] - v[1676] + v[1678] + v[2879] * v[723] + v[2878] * v[724] - 2e0*(v[1679] + v[2877])*v[725];
			v[2880] = (v[1676] + v[1678]) / 2e0;
			v[1663] = v[1663] + v[1677] - v[1681] + v[2880] * v[723] - 2e0*(v[1683] + v[2877])*v[724] + v[2878] * v[725];
			v[1662] = v[1662] - v[1680] + v[1682] - 2e0*(v[1679] + v[1683] - v[1684])*v[723] + v[2880] * v[724]
				+ v[2879] * v[725];
			v[2881] = v[1662] * v[711] + v[1663] * v[712] + v[1664] * v[713];
			v[1661] = v[1661] + v[2881] * v[716];
			v[1659] = v[1659] + v[2881] * v[722];
			v[1660] = v[1660] + v[1661] * v[721];
			v[1685] = v[1685] + 2e0*v[1659] * v[1699];
			v[1686] = v[1686] + (v[1685] * v[1700]) / 2e0;
			v[1658] = v[1658] + v[1660] + v[1686];
			v[2882] = v[1658] / v[714];
			v[1657] = v[1657] + v[1664] * v[2762] + v[2882] * v[713];
			v[1656] = v[1656] + v[1663] * v[2762] + v[2882] * v[712];
			v[1655] = v[1655] + v[1662] * v[2762] + v[2882] * v[711];
			v[1645] = v[1645] - v[1657] * v[381] + v[1656] * v[382];
			v[1653] = v[1653] - v[1656] * v[380] + v[1655] * v[381];
			v[1649] = v[1649] + v[1657] * v[380] - v[1655] * v[382];
		}
		else {
		};
		v[1969] = -(v[1864] * v[682]) - v[1863] * v[683] - v[1862] * v[684] + v[2885] * v[685] + v[2884] * v[687]
			+ v[2883] * v[689] - v[1857] * v[691] - v[1856] * v[692] - v[1855] * v[693];
		v[1970] = -(v[1864] * v[697]) - v[1863] * v[698] - v[1862] * v[699] + v[2885] * v[700] + v[2884] * v[702]
			+ v[2883] * v[704] - v[1857] * v[706] - v[1856] * v[707] - v[1855] * v[708];
		v[1971] = v[1864] * v[652] + v[1863] * v[653] + v[1862] * v[654] - v[2885] * v[655] - v[2884] * v[657] - v[2883] * v[659]
			+ v[1857] * v[661] + v[1856] * v[662] + v[1855] * v[663];
		v[1972] = v[1864] * v[667] + v[1863] * v[668] + v[1862] * v[669] - v[2885] * v[670] - v[2884] * v[672] - v[2883] * v[674]
			+ v[1857] * v[676] + v[1856] * v[677] + v[1855] * v[678];
		v[1973] = -(v[1851] * v[682]) - v[1850] * v[683] - v[1849] * v[684] + v[2888] * v[685] + v[2887] * v[687]
			+ v[2886] * v[689] - v[1844] * v[691] - v[1843] * v[692] - v[1842] * v[693];
		v[1974] = -(v[1851] * v[697]) - v[1850] * v[698] - v[1849] * v[699] + v[2888] * v[700] + v[2887] * v[702]
			+ v[2886] * v[704] - v[1844] * v[706] - v[1843] * v[707] - v[1842] * v[708];
		v[1975] = v[1851] * v[652] + v[1850] * v[653] + v[1849] * v[654] - v[2888] * v[655] - v[2887] * v[657] - v[2886] * v[659]
			+ v[1844] * v[661] + v[1843] * v[662] + v[1842] * v[663];
		v[1976] = v[1851] * v[667] + v[1850] * v[668] + v[1849] * v[669] - v[2888] * v[670] - v[2887] * v[672] - v[2886] * v[674]
			+ v[1844] * v[676] + v[1843] * v[677] + v[1842] * v[678];
		v[1977] = -(v[1838] * v[682]) - v[1837] * v[683] - v[1836] * v[684] + v[2891] * v[685] + v[2890] * v[687]
			+ v[2889] * v[689] - v[1831] * v[691] - v[1830] * v[692] - v[1829] * v[693];
		v[1978] = -(v[1838] * v[697]) - v[1837] * v[698] - v[1836] * v[699] + v[2891] * v[700] + v[2890] * v[702]
			+ v[2889] * v[704] - v[1831] * v[706] - v[1830] * v[707] - v[1829] * v[708];
		v[1979] = v[1838] * v[652] + v[1837] * v[653] + v[1836] * v[654] - v[2891] * v[655] - v[2890] * v[657] - v[2889] * v[659]
			+ v[1831] * v[661] + v[1830] * v[662] + v[1829] * v[663];
		v[1980] = v[1838] * v[667] + v[1837] * v[668] + v[1836] * v[669] - v[2891] * v[670] - v[2890] * v[672] - v[2889] * v[674]
			+ v[1831] * v[676] + v[1830] * v[677] + v[1829] * v[678];
		v[1984] = -((*ct)*(v[1138] * v[1797] + v[1139] * v[1798] + v[1140] * v[1799])) + (*cn)*(v[1794] * v[1981]
			+ v[1795] * v[1982] + v[1796] * v[1983]) - v[1788] * v[663] - v[1787] * v[678] + v[1790] * v[693] + v[1789] * v[708];
		v[2086] = v[1984] * v[2708];
		v[1988] = -((*ct)*(v[1134] * v[1797] + v[1135] * v[1798] + v[1136] * v[1799])) + (*cn)*(v[1794] * v[1985]
			+ v[1795] * v[1986] + v[1796] * v[1987]) - v[1788] * v[662] - v[1787] * v[677] + v[1790] * v[692] + v[1789] * v[707];
		v[2898] = v[1988] * v[333];
		v[2090] = v[1988] * v[2707];
		v[1992] = -((*ct)*(v[1130] * v[1797] + v[1131] * v[1798] + v[1132] * v[1799])) + (*cn)*(v[1794] * v[1989]
			+ v[1795] * v[1990] + v[1796] * v[1991]) - v[1788] * v[661] - v[1787] * v[676] + v[1790] * v[691] + v[1789] * v[706];
		v[2901] = v[1992] * v[333];
		v[2093] = v[1992] * v[2706];
		v[1996] = -((*ct)*(v[1114] * v[1797] + v[1115] * v[1798] + v[1116] * v[1799])) + (*cn)*(v[1794] * v[1993]
			+ v[1795] * v[1994] + v[1796] * v[1995]) - v[1788] * v[654] - v[1787] * v[669] + v[1790] * v[684] + v[1789] * v[699];
		v[2155] = v[1996] * v[2705];
		v[2000] = -((*ct)*(v[1110] * v[1797] + v[1111] * v[1798] + v[1112] * v[1799])) + (*cn)*(v[1794] * v[1997]
			+ v[1795] * v[1998] + v[1796] * v[1999]) - v[1788] * v[653] - v[1787] * v[668] + v[1790] * v[683] + v[1789] * v[698];
		v[2906] = v[2000] * v[307];
		v[2159] = v[2000] * v[2704];
		v[2004] = -((*ct)*(v[1106] * v[1797] + v[1107] * v[1798] + v[1108] * v[1799])) + (*cn)*(v[1794] * v[2001]
			+ v[1795] * v[2002] + v[1796] * v[2003]) - v[1788] * v[652] - v[1787] * v[667] + v[1790] * v[682] + v[1789] * v[697];
		v[5088] = v[2965];
		v[5089] = -v[2966];
		v[5090] = v[2967];
		v[5091] = v[1170] * v[1580] + v[1169] * v[1584] + v[1497] * (v[2165] + v[1169] * v[2259] + v[1170] * v[2261])
			+ v[2004] * v[2293] + v[2000] * v[2710] + v[1996] * v[2712] + v[1486] * v[2968] + v[1489] * v[2969] + v[1492] * v[2970]
			+ v[225] * v[2971] + v[230] * v[2972] + v[221] * v[2973];
		v[5092] = v[1171] * v[1580] + v[1169] * v[1586] + v[1417] * v[1996] + v[1418] * v[2004] + v[1497] * (v[2158]
			+ v[1169] * v[2258] + v[1171] * v[2262]) + v[2000] * v[2294] + v[1493] * v[2806] + v[1487] * v[2808] + v[1490] * v[2974]
			+ v[222] * v[2975] + v[227] * v[2976] + v[1594] * v[3004];
		v[5093] = v[1171] * v[1584] + v[1170] * v[1586] + v[1419] * v[2000] + v[1420] * v[2004] + v[1497] * (v[2150]
			+ v[1170] * v[2260] + v[1171] * v[2263]) + v[1996] * v[2295] + v[1491] * v[2807] + v[1488] * v[2977] + v[1494] * v[2978]
			+ v[233] * v[2979] + v[1598] * v[3002] + v[1602] * v[3003];
		v[5094] = -v[2965];
		v[5095] = v[2966];
		v[5096] = -v[2967];
		v[5097] = v[1167] * v[1608] + v[1166] * v[1612] + v[1548] * (v[2096] + v[1166] * v[2230] + v[1167] * v[2232])
			+ v[1992] * v[2296] + v[1988] * v[2715] + v[1984] * v[2717] + v[1537] * v[2980] + v[1540] * v[2981] + v[1543] * v[2982]
			+ v[265] * v[2983] + v[270] * v[2984] + v[261] * v[2985];
		v[5098] = v[1168] * v[1608] + v[1166] * v[1614] + v[1423] * v[1984] + v[1424] * v[1992] + v[1548] * (v[2089]
			+ v[1166] * v[2229] + v[1168] * v[2233]) + v[1988] * v[2297] + v[1544] * v[2793] + v[1538] * v[2795] + v[1541] * v[2986]
			+ v[262] * v[2987] + v[267] * v[2988] + v[1622] * v[3001];
		v[5099] = v[1168] * v[1612] + v[1167] * v[1614] + v[1425] * v[1988] + v[1426] * v[1992] + v[1548] * (v[2081]
			+ v[1167] * v[2231] + v[1168] * v[2234]) + v[1984] * v[2298] + v[1542] * v[2794] + v[1539] * v[2989] + v[1545] * v[2990]
			+ v[273] * v[2991] + v[1626] * v[2999] + v[1630] * v[3000];
		v[2909] = v[2004] * v[307];
		v[2162] = v[2004] * v[2703];
		v[1653] = v[1653] + 2e0*v[1871] * v[373];
		v[1649] = v[1649] + 2e0*v[1872] * v[372];
		v[1645] = v[1645] + 2e0*v[1881] * v[371];
		v[2893] = (-(v[1263] * v[2005] * v[2892]) + v[1642] * v[2006] * v[368] + (v[1179] * v[1634] + v[1180] * v[1637]
			+ v[1181] * v[1640] + v[1645] * v[356] + v[1649] * v[357] + v[1653] * v[358])*v[369]) / v[1264];
		v[2009] = (*epsn)*v[1796] + v[1640] * v[2792] + v[1181] * v[2827] + v[2893] * v[358] + v[1653] * v[370];
		v[2011] = (*epsn)*v[1795] + v[1637] * v[2792] + v[1180] * v[2827] + v[2893] * v[357] + v[1649] * v[370];
		v[2013] = (*epsn)*v[1794] + v[1634] * v[2792] + v[1179] * v[2827] + v[2893] * v[356] + v[1645] * v[370];
		v[1953] = v[1953] - v[2009];
		v[1957] = v[1957] + v[2009];
		v[1945] = v[1945] - v[2011];
		v[1949] = v[1949] + v[2011];
		v[1937] = v[1937] - v[2013];
		v[1941] = v[1941] + v[2013];
		v[2014] = v[1166] * v[2894] + v[1984] * v[342];
		v[2015] = v[1166] * v[2895] + v[1984] * v[336];
		v[2016] = v[1283] * v[1984] + v[1166] * (v[1548] * v[2899] + v[1565] / v[333]) + (v[2014] * v[267]) / v[333];
		v[2017] = v[1290] * v[1984] + v[1166] * (v[1548] * v[3139] + v[1557] / v[333]) + (v[2015] * v[261]) / v[333];
		v[2018] = (v[2015] * v[270] + v[2014] * v[272] + v[1166] * (v[1573] + v[2826] * v[3140]) + v[1279] * v[1984] * v[333])
			/ v[333];
		v[2019] = v[1167] * v[2896] + v[1988] * v[344];
		v[2020] = v[1167] * v[2895] + v[1988] * v[336];
		v[2023] = v[1278] * v[1988] + v[1167] * (v[1548] * v[2897] + v[1571] / v[333]) + (v[2019] * v[273]) / v[333];
		v[2024] = v[2014] + v[2019];
		v[2025] = v[2016] + v[2023];
		v[2027] = (v[1281] * v[1539] + v[2020] * v[261] + v[2024] * v[263] + v[2225] * v[2826] + v[1289] * v[2898] + v[1167] *
			(v[1553] + v[2826] * v[2900])) / v[333];
		v[2030] = (v[2020] * v[265] + v[2019] * v[268] + v[1282] * v[2898] + v[1167] * (v[1562] + v[2826] * v[3141])) / v[333];
		v[2031] = v[1168] * v[2894] + v[1992] * v[342];
		v[2032] = v[1168] * v[2896] + v[1992] * v[344];
		v[2033] = v[2020] + v[2031];
		v[2036] = (v[2031] * v[262] + v[2032] * v[263] + v[1287] * v[2901] + v[1168] * (v[1554] + v[2826] * v[3142])) / v[333];
		v[2037] = v[2015] + v[2032];
		v[2039] = (v[1292] * v[1542] + v[2031] * v[267] + v[2037] * v[268] + v[2224] * v[2826] + v[1294] * v[2901] + v[1168] *
			(v[1561] + v[2826] * v[3143])) / v[333];
		v[2040] = v[2027] + v[2039];
		v[2042] = (v[1293] * v[1544] + v[2033] * v[272] + v[2032] * v[273] + v[2223] * v[2826] + v[1298] * v[2901] + v[1168] *
			(v[1569] + v[2826] * v[3144])) / v[333];
		v[2043] = v[2017] + v[2042];
		v[2044] = v[1169] * v[2902] + v[1996] * v[316];
		v[2045] = v[1169] * v[2903] + v[1996] * v[310];
		v[2046] = v[1310] * v[1996] + v[1169] * (v[1497] * v[2907] + v[1514] / v[307]) + (v[2044] * v[227]) / v[307];
		v[2047] = v[1317] * v[1996] + (v[2045] * v[221]) / v[307] + v[1169] * (v[1506] / v[307] + v[1497] * v[3145]);
		v[2048] = (v[2045] * v[230] + v[2044] * v[232] + v[1306] * v[1996] * v[307] + v[1169] * (v[1522] + v[2825] * v[3146]))
			/ v[307];
		v[2049] = v[1170] * v[2904] + v[2000] * v[318];
		v[2050] = v[1170] * v[2903] + v[2000] * v[310];
		v[2053] = v[1305] * v[2000] + v[1170] * (v[1497] * v[2905] + v[1520] / v[307]) + (v[2049] * v[233]) / v[307];
		v[2054] = v[2044] + v[2049];
		v[2055] = v[2046] + v[2053];
		v[2057] = (v[1308] * v[1488] + v[2050] * v[221] + v[2054] * v[223] + v[2254] * v[2825] + v[1316] * v[2906] + v[1170] *
			(v[1502] + v[2825] * v[2908])) / v[307];
		v[2060] = (v[2050] * v[225] + v[2049] * v[228] + v[1309] * v[2906] + v[1170] * (v[1511] + v[2825] * v[3147])) / v[307];
		v[2061] = v[1171] * v[2902] + v[2004] * v[316];
		v[2062] = v[1171] * v[2904] + v[2004] * v[318];
		v[2063] = v[2050] + v[2061];
		v[2066] = (v[2061] * v[222] + v[2062] * v[223] + v[1314] * v[2909] + v[1171] * (v[1503] + v[2825] * v[3148])) / v[307];
		v[2067] = v[2045] + v[2062];
		v[2069] = (v[1319] * v[1491] + v[2061] * v[227] + v[2067] * v[228] + v[2253] * v[2825] + v[1321] * v[2909] + v[1171] *
			(v[1510] + v[2825] * v[3149])) / v[307];
		v[2070] = v[2057] + v[2069];
		v[2072] = (v[1320] * v[1493] + v[2063] * v[232] + v[2062] * v[233] + v[2252] * v[2825] + v[1325] * v[2909] + v[1171] *
			(v[1518] + v[2825] * v[3150])) / v[307];
		v[2073] = v[2047] + v[2072];
		v[1954] = v[1954] + dGBp[2][0] * v[1977] + dGBp[2][1] * v[1978] - GBp[2] * v[2009];
		v[1955] = v[1955] + dGBp[1][0] * v[1977] + dGBp[1][1] * v[1978] - GBp[1] * v[2009];
		v[1956] = v[1956] + dGBp[0][0] * v[1977] + dGBp[0][1] * v[1978] - GBp[0] * v[2009];
		v[1946] = v[1946] + dGBp[2][0] * v[1973] + dGBp[2][1] * v[1974] - GBp[2] * v[2011];
		v[1947] = v[1947] + dGBp[1][0] * v[1973] + dGBp[1][1] * v[1974] - GBp[1] * v[2011];
		v[1948] = v[1948] + dGBp[0][0] * v[1973] + dGBp[0][1] * v[1974] - GBp[0] * v[2011];
		v[2074] = v[1172] * v[1560] + v[1159] * v[1568] + v[1154] * v[1576] + v[1969] * v[279] + v[1973] * v[282]
			+ v[1977] * v[285];
		v[2075] = v[1172] * v[1559] + v[1159] * v[1567] + v[1154] * v[1575] + v[1969] * v[278] + v[1973] * v[281]
			+ v[1977] * v[284];
		v[2076] = v[1172] * v[1558] + v[1159] * v[1566] + v[1154] * v[1574] + v[1969] * v[277] + v[1973] * v[280]
			+ v[1977] * v[283];
		v[1938] = v[1938] + dGBp[2][0] * v[1969] + dGBp[2][1] * v[1970] - GBp[2] * v[2013];
		v[1939] = v[1939] + dGBp[1][0] * v[1969] + dGBp[1][1] * v[1970] - GBp[1] * v[2013];
		v[1940] = v[1940] + dGBp[0][0] * v[1969] + dGBp[0][1] * v[1970] - GBp[0] * v[2013];
		v[2077] = v[1173] * v[1560] + v[1161] * v[1568] + v[1155] * v[1576] + v[1970] * v[279] + v[1974] * v[282]
			+ v[1978] * v[285];
		v[2078] = v[1173] * v[1559] + v[1161] * v[1567] + v[1155] * v[1575] + v[1970] * v[278] + v[1974] * v[281]
			+ v[1978] * v[284];
		v[2079] = v[1173] * v[1558] + v[1161] * v[1566] + v[1155] * v[1574] + v[1970] * v[277] + v[1974] * v[280]
			+ v[1978] * v[283];
		v[2080] = QBi[2][2] * v[1954] + QBi[2][1] * v[1955] + QBi[2][0] * v[1956] + v[1423] * v[2019] + v[2032] * v[2717] +
			(v[2991] + v[1548] * v[3151])*v[344] + v[2086] * v[346];
		v[2083] = QBi[1][2] * v[1954] + QBi[1][1] * v[1955] + QBi[1][0] * v[1956] + v[2014] * v[2298] + v[1622] * v[2707]
			+ v[2033] * v[2717] + v[1548] * v[3152] + (v[1293] * v[1605]) / v[333] + v[2090] * v[343];
		v[2085] = QBi[0][2] * v[1954] + QBi[0][1] * v[1955] + QBi[0][0] * v[1956] + v[2015] * v[2298] + v[2093] * v[331] +
			(v[2984] + v[1548] * v[3153])*v[336];
		v[2087] = QBi[2][2] * v[1946] + QBi[2][1] * v[1947] + QBi[2][0] * v[1948] + v[2019] * v[2297] + v[1630] * v[2708]
			+ v[2037] * v[2715] + v[2086] * v[2718] + v[1548] * v[3154] + (v[1292] * v[1609]) / v[333];
		v[2091] = QBi[1][2] * v[1946] + QBi[1][1] * v[1947] + QBi[1][0] * v[1948] + v[1425] * v[2014] + v[2031] * v[2715]
			+ v[2090] * v[337] + (v[2988] + v[1548] * v[3155])*v[342];
		v[2094] = QBi[0][2] * v[1946] + QBi[0][1] * v[1947] + QBi[0][0] * v[1948] + v[2020] * v[2297] + v[2093] * v[332] +
			(v[2983] + v[1548] * v[3156])*v[336];
		v[2095] = QBi[2][2] * v[1938] + QBi[2][1] * v[1939] + QBi[2][0] * v[1940] + v[1424] * v[2024] + v[2032] * v[2296]
			+ v[1626] * v[2708] + v[2086] * v[2716] + v[1548] * v[3157] + (v[1281] * v[1610]) / v[333];
		v[2097] = QBi[1][2] * v[1938] + QBi[1][1] * v[1939] + QBi[1][0] * v[1940] + v[2031] * v[2296] + v[2090] * v[2714] +
			(v[1548] * v[2912] + v[2987])*v[342];
		v[2100] = QBi[0][2] * v[1938] + QBi[0][1] * v[1939] + QBi[0][0] * v[1940] + v[1426] * v[2015] + v[1424] * v[2020]
			+ v[2093] * v[334] + (v[2985] + v[1548] * v[3158])*v[336];
		v[2101] = -(v[1899] * v[856]);
		v[2102] = -(v[1899] * v[854]);
		v[2103] = -(v[1899] * v[852]);
		v[2104] = -(v[1902] * v[856]);
		v[2105] = -(v[1902] * v[852]);
		v[2106] = -(v[1902] * v[854]);
		v[2107] = -(v[1904] * v[854]);
		v[2108] = -(v[1904] * v[852]);
		v[2109] = -(v[1904] * v[856]);
		v[2110] = -(v[1905] * v[856]);
		v[2111] = -(v[1905] * v[854]);
		v[2112] = -(v[1905] * v[852]);
		v[2113] = -(v[1908] * v[852]);
		v[2114] = -(v[1908] * v[856]);
		v[2115] = -(v[1908] * v[854]);
		v[2116] = -(v[1910] * v[854]);
		v[2117] = -(v[1910] * v[856]);
		v[2118] = -(v[1910] * v[852]);
		v[2119] = -(v[1285] * v[1527]) + 2e0*v[1284] * v[1529] - v[1296] * v[1530] + v[2107] + v[2110] + v[2113] + v[2116]
			+ 2e0*v[2030] * v[517] - v[2040] * v[520] - v[2025] * v[523];
		v[2120] = -(v[1911] * v[856]);
		v[2121] = -(v[1911] * v[852]);
		v[2122] = -(v[1911] * v[854]);
		v[2123] = v[1300] * v[1527] - v[1296] * v[1529] + 2e0*v[1291] * v[1530] - v[2102] - v[2105] - v[2117] - v[2120]
			- v[2040] * v[517] + 2e0*v[2036] * v[520] + v[2043] * v[523];
		v[2124] = -2e0*v[1349] * v[1528] + v[2097] * v[258] - v[2107] * v[510] + v[2102] * v[511] - v[2106] * v[512];
		v[2125] = -(v[1913] * v[856]);
		v[2126] = -(v[1913] * v[854]);
		v[2127] = -(v[1913] * v[852]);
		v[2128] = -(v[1914] * v[854]);
		v[2129] = -(v[1914] * v[856]);
		v[2130] = -(v[1914] * v[852]);
		v[2131] = -(v[1268] * v[1560]) - v[1267] * v[1568] - v[1265] * v[1576] - v[2013] * v[279] - v[2011] * v[282]
			- v[2009] * v[285] + v[1904] * v[561] + v[1899] * v[562] + v[1902] * v[563] + v[1905] * v[576] + v[1911] * v[577]
			+ v[1908] * v[578] + v[1913] * v[591] + v[1914] * v[592] + v[1910] * v[593];
		v[2132] = -(v[1268] * v[1559]) - v[1267] * v[1567] - v[1265] * v[1575] - v[2013] * v[278] - v[2011] * v[281]
			- v[2009] * v[284] + v[1904] * v[558] + v[1899] * v[559] + v[1902] * v[560] + v[1905] * v[573] + v[1911] * v[574]
			+ v[1908] * v[575] + v[1913] * v[588] + v[1914] * v[589] + v[1910] * v[590];
		v[2133] = -(v[1268] * v[1558]) - v[1267] * v[1566] - v[1265] * v[1574] - v[2013] * v[277] - v[2011] * v[280]
			- v[2009] * v[283] + v[1904] * v[555] + v[1899] * v[556] + v[1902] * v[557] + v[1905] * v[570] + v[1911] * v[571]
			+ v[1908] * v[572] + v[1913] * v[585] + v[1914] * v[586] + v[1910] * v[587];
		v[2134] = 2e0*v[1277] * v[1527] - v[1285] * v[1529] + v[1300] * v[1530] - v[2108] - v[2121] - v[2125] - v[2128]
			- v[2025] * v[517] + v[2043] * v[520] + 2e0*v[2018] * v[523];
		v[2135] = -2e0*v[1347] * v[1528] + v[2095] * v[258] - v[2108] * v[510] + v[2103] * v[511] - v[2105] * v[512];
		v[2136] = -2e0*v[1340] * v[1528] + v[2094] * v[258] - v[2110] * v[510] + v[2120] * v[511] - v[2114] * v[512];
		v[2137] = v[2104] + v[2115];
		v[2138] = -2e0*v[1335] * v[1528] + v[2087] * v[258] - v[2112] * v[510] + v[2121] * v[511] - v[2113] * v[512];
		v[2139] = -2e0*v[1334] * v[1528] + v[2085] * v[258] - v[2125] * v[510] + v[2129] * v[511] - v[2117] * v[512];
		v[2140] = -2e0*v[1332] * v[1528] + v[2083] * v[258] - v[2126] * v[510] + v[2128] * v[511] - v[2116] * v[512];
		v[2141] = v[2111] + v[2127];
		v[2142] = v[2101] + v[2130];
		v[1958] = v[1958] + dGAp[2][0] * v[1979] + dGAp[2][1] * v[1980] + GAp[2] * v[2009];
		v[1959] = v[1959] + dGAp[1][0] * v[1979] + dGAp[1][1] * v[1980] + GAp[1] * v[2009];
		v[1960] = v[1960] + dGAp[0][0] * v[1979] + dGAp[0][1] * v[1980] + GAp[0] * v[2009];
		v[1950] = v[1950] + dGAp[2][0] * v[1975] + dGAp[2][1] * v[1976] + GAp[2] * v[2011];
		v[1951] = v[1951] + dGAp[1][0] * v[1975] + dGAp[1][1] * v[1976] + GAp[1] * v[2011];
		v[1952] = v[1952] + dGAp[0][0] * v[1975] + dGAp[0][1] * v[1976] + GAp[0] * v[2011];
		v[2143] = v[1174] * v[1507] + v[1163] * v[1515] + v[1156] * v[1523] + v[1971] * v[237] + v[1975] * v[240]
			+ v[1979] * v[243];
		v[2144] = v[1174] * v[1508] + v[1163] * v[1516] + v[1156] * v[1524] + v[1971] * v[238] + v[1975] * v[241]
			+ v[1979] * v[244];
		v[2145] = v[1174] * v[1509] + v[1163] * v[1517] + v[1156] * v[1525] + v[1971] * v[239] + v[1975] * v[242]
			+ v[1979] * v[245];
		v[1942] = v[1942] + dGAp[2][0] * v[1971] + dGAp[2][1] * v[1972] + GAp[2] * v[2013];
		v[1943] = v[1943] + dGAp[1][0] * v[1971] + dGAp[1][1] * v[1972] + GAp[1] * v[2013];
		v[1944] = v[1944] + dGAp[0][0] * v[1971] + dGAp[0][1] * v[1972] + GAp[0] * v[2013];
		v[2146] = v[1175] * v[1507] + v[1165] * v[1515] + v[1157] * v[1523] + v[1972] * v[237] + v[1976] * v[240]
			+ v[1980] * v[243];
		v[2147] = v[1175] * v[1508] + v[1165] * v[1516] + v[1157] * v[1524] + v[1972] * v[238] + v[1976] * v[241]
			+ v[1980] * v[244];
		v[2148] = v[1175] * v[1509] + v[1165] * v[1517] + v[1157] * v[1525] + v[1972] * v[239] + v[1976] * v[242]
			+ v[1980] * v[245];
		v[2149] = QAi[2][2] * v[1958] + QAi[2][1] * v[1959] + QAi[2][0] * v[1960] + v[1417] * v[2049] + v[2062] * v[2712] +
			(v[2979] + v[1497] * v[3159])*v[318] + v[2155] * v[320];
		v[2152] = QAi[1][2] * v[1958] + QAi[1][1] * v[1959] + QAi[1][0] * v[1960] + v[2044] * v[2295] + v[1594] * v[2704]
			+ v[2063] * v[2712] + (v[1320] * v[1577]) / v[307] + v[1497] * v[3160] + v[2159] * v[317];
		v[2154] = QAi[0][2] * v[1958] + QAi[0][1] * v[1959] + QAi[0][0] * v[1960] + v[2045] * v[2295] + v[2162] * v[305]
			+ v[310] * (v[2972] + v[1497] * v[3161]);
		v[2156] = QAi[2][2] * v[1950] + QAi[2][1] * v[1951] + QAi[2][0] * v[1952] + v[2049] * v[2294] + v[1602] * v[2705]
			+ v[2067] * v[2710] + v[2155] * v[2713] + (v[1319] * v[1581]) / v[307] + v[1497] * v[3162];
		v[2160] = QAi[1][2] * v[1950] + QAi[1][1] * v[1951] + QAi[1][0] * v[1952] + v[1419] * v[2044] + v[2061] * v[2710]
			+ v[2159] * v[311] + v[316] * (v[2976] + v[1497] * v[3163]);
		v[2163] = QAi[0][2] * v[1950] + QAi[0][1] * v[1951] + QAi[0][0] * v[1952] + v[2050] * v[2294] + v[2162] * v[306]
			+ v[310] * (v[2971] + v[1497] * v[3164]);
		v[2164] = QAi[2][2] * v[1942] + QAi[2][1] * v[1943] + QAi[2][0] * v[1944] + v[1418] * v[2054] + v[2062] * v[2293]
			+ v[1598] * v[2705] + v[2155] * v[2711] + (v[1308] * v[1582]) / v[307] + v[1497] * v[3165];
		v[2166] = QAi[1][2] * v[1942] + QAi[1][1] * v[1943] + QAi[1][0] * v[1944] + v[2061] * v[2293] + v[2159] * v[2709] +
			(v[1497] * v[2915] + v[2975])*v[316];
		v[2169] = QAi[0][2] * v[1942] + QAi[0][1] * v[1943] + QAi[0][0] * v[1944] + v[1420] * v[2045] + v[1418] * v[2050]
			+ v[2162] * v[308] + v[310] * (v[2973] + v[1497] * v[3166]);
		v[2170] = v[1916] * v[868];
		v[2171] = v[1916] * v[866];
		v[2172] = v[1916] * v[864];
		v[2173] = v[1918] * v[868];
		v[2174] = v[1918] * v[864];
		v[2175] = v[1918] * v[866];
		v[2176] = v[1920] * v[866];
		v[2177] = v[1920] * v[864];
		v[2178] = v[1920] * v[868];
		v[2179] = v[1921] * v[868];
		v[2180] = v[1921] * v[866];
		v[2181] = v[1921] * v[864];
		v[2182] = v[1922] * v[864];
		v[2183] = v[1922] * v[868];
		v[2184] = v[1922] * v[866];
		v[2185] = v[1923] * v[866];
		v[2186] = v[1923] * v[868];
		v[2187] = v[1923] * v[864];
		v[2188] = -(v[1312] * v[1476]) + 2e0*v[1311] * v[1478] - v[1323] * v[1479] + v[2176] + v[2179] + v[2182] + v[2185]
			+ 2e0*v[2060] * v[418] - v[2070] * v[421] - v[2055] * v[424];
		v[2189] = v[1924] * v[868];
		v[2190] = v[1924] * v[864];
		v[2191] = v[1924] * v[866];
		v[2192] = v[1327] * v[1476] - v[1323] * v[1478] + 2e0*v[1318] * v[1479] - v[2171] - v[2174] - v[2186] - v[2189]
			- v[2070] * v[418] + 2e0*v[2066] * v[421] + v[2073] * v[424];
		v[2193] = -2e0*v[1387] * v[1477] + v[2166] * v[218] - v[2176] * v[411] + v[2171] * v[412] - v[2175] * v[413];
		v[2194] = v[1925] * v[868];
		v[2195] = v[1925] * v[866];
		v[2196] = v[1925] * v[864];
		v[2197] = v[1926] * v[866];
		v[2198] = v[1926] * v[868];
		v[2199] = v[1926] * v[864];
		v[2200] = v[1268] * v[1507] + v[1267] * v[1515] + v[1265] * v[1523] + v[2013] * v[237] + v[2011] * v[240]
			+ v[2009] * v[243] + v[1920] * v[456] + v[1916] * v[457] + v[1918] * v[458] + v[1921] * v[471] + v[1924] * v[472]
			+ v[1922] * v[473] + v[1925] * v[486] + v[1926] * v[487] + v[1923] * v[488];
		v[2201] = v[1268] * v[1508] + v[1267] * v[1516] + v[1265] * v[1524] + v[2013] * v[238] + v[2011] * v[241]
			+ v[2009] * v[244] + v[1920] * v[459] + v[1916] * v[460] + v[1918] * v[461] + v[1921] * v[474] + v[1924] * v[475]
			+ v[1922] * v[476] + v[1925] * v[489] + v[1926] * v[490] + v[1923] * v[491];
		v[2202] = v[1268] * v[1509] + v[1267] * v[1517] + v[1265] * v[1525] + v[2013] * v[239] + v[2011] * v[242]
			+ v[2009] * v[245] + v[1920] * v[462] + v[1916] * v[463] + v[1918] * v[464] + v[1921] * v[477] + v[1924] * v[478]
			+ v[1922] * v[479] + v[1925] * v[492] + v[1926] * v[493] + v[1923] * v[494];
		v[2203] = 2e0*v[1304] * v[1476] - v[1312] * v[1478] + v[1327] * v[1479] - v[2177] - v[2190] - v[2194] - v[2197]
			- v[2055] * v[418] + v[2073] * v[421] + 2e0*v[2048] * v[424];
		v[2204] = -2e0*v[1385] * v[1477] + v[2164] * v[218] - v[2177] * v[411] + v[2172] * v[412] - v[2174] * v[413];
		v[2205] = -2e0*v[1378] * v[1477] + v[2163] * v[218] - v[2179] * v[411] + v[2189] * v[412] - v[2183] * v[413];
		v[2206] = v[2173] + v[2184];
		v[2207] = -2e0*v[1373] * v[1477] + v[2156] * v[218] - v[2181] * v[411] + v[2190] * v[412] - v[2182] * v[413];
		v[2208] = -2e0*v[1372] * v[1477] + v[2154] * v[218] - v[2194] * v[411] + v[2198] * v[412] - v[2186] * v[413];
		v[2209] = -2e0*v[1370] * v[1477] + v[2152] * v[218] - v[2195] * v[411] + v[2197] * v[412] - v[2185] * v[413];
		v[2210] = v[2180] + v[2196];
		v[2211] = v[2170] + v[2199];
		v[2212] = (-2e0*v[2016] + 2e0*v[2023] - v[2109] * v[513] - 2e0*v[2107] * v[518] - 2e0*v[2108] * v[524]
			- 2e0*v[2110] * v[528] - v[2111] * v[532] - 2e0*v[2112] * v[537] - 2e0*v[2125] * v[541] - 2e0*v[2126] * v[545]
			- v[2127] * v[550]) / 2e0;
		v[2213] = (-2e0*v[1352] * v[1528] + v[2100] * v[258] - v[2109] * v[510] + v[2101] * v[511] - v[2104] * v[512]) / 2e0;
		v[2215] = -v[2017] + v[2042] + (v[2101] * v[513]) / 2e0 + v[2102] * v[518] + v[2103] * v[524] + v[2120] * v[528] +
			(v[2122] * v[532]) / 2e0 + v[2121] * v[537] + v[2129] * v[541] + v[2128] * v[545] + (v[2130] * v[550]) / 2e0;
		v[2216] = (-2e0*v[1338] * v[1528] + v[2091] * v[258] - v[2111] * v[510] + v[2122] * v[511] - v[2115] * v[512]) / 2e0;
		v[2217] = (-2e0*v[2027] + 2e0*v[2039] - v[2104] * v[513] - 2e0*v[2106] * v[518] - 2e0*v[2105] * v[524]
			- 2e0*v[2114] * v[528] - v[2115] * v[532] - 2e0*v[2113] * v[537] - 2e0*v[2117] * v[541] - 2e0*v[2116] * v[545]
			- v[2118] * v[550]) / 2e0;
		v[2282] = v[1531] * v[2215] - v[2216] + v[1533] * v[2217] + 24e0*v[1448] * v[3167] * v[3168] - 2e0*v[1413] * (
			-4e0*v[1366] * v[1448] + 4e0*v[2212] * v[405] + v[4979 + i1143]) - 2e0*(v[1352] * v[1449] + v[1329] * v[1450]
				+ 2e0*v[1349] * v[1469] + 2e0*v[1340] * v[1470] + 2e0*v[1335] * v[1471] + 2e0*v[1332] * v[1472]
				+ 2e0*v[1347] * v[1473] + 2e0*v[1334] * v[1474] + v[1338] * v[1526] + 2e0*v[2103] - 2e0*v[2106] - 2e0*v[2112]
				+ 2e0*v[2114] + alphaB[1] * v[2119] - alphaB[0] * v[2123] + 2e0*v[2126] - 2e0*v[2129] - alphaB[2] * v[2134]
				- 8e0*v[1528] * v[2219] + 4e0*v[258] * (v[2018] + v[1622] * v[2021] + v[1619] * v[2022] + v[1613] * v[2024]
					+ v[1628] * v[2026] + v[1632] * v[2028] + v[1630] * v[2029] + v[2030] + v[1607] * v[2033] + v[1615] * v[2034]
					+ v[1626] * v[2035] + v[2036] + v[1611] * v[2037] + v[1618] * v[2038] + v[1624] * v[2041] + v[1573] * v[2081]
					+ v[1571] * v[2082] + v[1569] * v[2084] + v[1565] * v[2088] + v[1562] * v[2089] + v[1561] * v[2092] + v[1554] * v[2096]
					+ v[1557] * v[2098] + v[1553] * v[2099] + v[1539] * v[2220] + v[1542] * v[2221] + v[1544] * v[2222] + v[1605] * v[2223]
					+ v[1609] * v[2224] + v[1610] * v[2225] + v[2014] * v[2229] + v[2015] * v[2230] + v[2019] * v[2231] + v[2020] * v[2232]
					+ v[2031] * v[2233] + v[2032] * v[2234] + v[1984] * v[2926] + v[1988] * v[2927] + v[1992] * v[2928]
					- 2e0*v[1548] * v[3171] * v[3172]) - v[2141] * v[405] - v[2142] * v[408] - v[2137] * v[410] + 2e0*v[4955 + i1143]
				+ v[2100] * v[513] + 2e0*v[2097] * v[518] + 2e0*v[2095] * v[524] + 2e0*v[2094] * v[528] + v[2091] * v[532]
				+ 2e0*v[2087] * v[537] + 2e0*v[2085] * v[541] + 2e0*v[2083] * v[545] + v[2080] * v[550])*v[977];
		v[2964] = v[2282] + (2e0*v[1329] * v[1528] - v[2080] * v[258] + v[2127] * v[510] - v[2130] * v[511] + v[2118] * v[512])
			/ 2e0;
		v[2962] = (v[2135] + v[2139]) / 2e0;
		v[2237] = v[2138] + v[2140];
		v[2238] = v[2124] + v[2136];
		v[2239] = ddGBp[2][0][1] * v[2074] + ddGBp[1][0][1] * v[2075] + ddGBp[0][0][1] * v[2076]
			+ ddGBp[2][1][1] * v[2077] + ddGBp[1][1][1] * v[2078] + ddGBp[0][1][1] * v[2079] + dGBp[2][1] * v[2131]
			+ dGBp[1][1] * v[2132] + dGBp[0][1] * v[2133];
		v[2240] = ddGBp[2][0][0] * v[2074] + ddGBp[1][0][0] * v[2075] + ddGBp[0][0][0] * v[2076]
			+ ddGBp[2][1][0] * v[2077] + ddGBp[1][1][0] * v[2078] + ddGBp[0][1][0] * v[2079] + dGBp[2][0] * v[2131]
			+ dGBp[1][0] * v[2132] + dGBp[0][0] * v[2133];
		v[2241] = (-2e0*v[2046] + 2e0*v[2053] - v[2178] * v[414] - 2e0*v[2176] * v[419] - 2e0*v[2177] * v[425]
			- 2e0*v[2179] * v[429] - v[2180] * v[433] - 2e0*v[2181] * v[438] - 2e0*v[2194] * v[442] - 2e0*v[2195] * v[446]
			- v[2196] * v[451]) / 2e0;
		v[2242] = (-2e0*v[1390] * v[1477] + v[2169] * v[218] - v[2178] * v[411] + v[2170] * v[412] - v[2173] * v[413]) / 2e0;
		v[2244] = -v[2047] + v[2072] + (v[2170] * v[414]) / 2e0 + v[2171] * v[419] + v[2172] * v[425] + v[2189] * v[429] +
			(v[2191] * v[433]) / 2e0 + v[2190] * v[438] + v[2198] * v[442] + v[2197] * v[446] + (v[2199] * v[451]) / 2e0;
		v[2245] = (-2e0*v[1376] * v[1477] + v[2160] * v[218] - v[2180] * v[411] + v[2191] * v[412] - v[2184] * v[413]) / 2e0;
		v[2246] = (-2e0*v[2057] + 2e0*v[2069] - v[2173] * v[414] - 2e0*v[2175] * v[419] - 2e0*v[2174] * v[425]
			- 2e0*v[2183] * v[429] - v[2184] * v[433] - 2e0*v[2182] * v[438] - 2e0*v[2186] * v[442] - 2e0*v[2185] * v[446]
			- v[2187] * v[451]) / 2e0;
		v[2274] = v[1480] * v[2244] - v[2245] + v[1482] * v[2246] + 24e0*v[1444] * v[3173] * v[3174] - 2e0*v[1411] * (
			-4e0*v[1404] * v[1444] + 4e0*v[2241] * v[399] + v[4991 + i1143]) - 2e0*(v[1390] * v[1445] + v[1367] * v[1446]
				+ 2e0*v[1387] * v[1463] + 2e0*v[1378] * v[1464] + 2e0*v[1373] * v[1465] + 2e0*v[1370] * v[1466]
				+ 2e0*v[1385] * v[1467] + 2e0*v[1372] * v[1468] + v[1376] * v[1475] + 2e0*v[2172] - 2e0*v[2175] - 2e0*v[2181]
				+ 2e0*v[2183] + alphaA[1] * v[2188] - alphaA[0] * v[2192] + 2e0*v[2195] - 2e0*v[2198] - alphaA[2] * v[2203]
				- 8e0*v[1477] * v[2248] + 4e0*v[218] * (v[2048] + v[1594] * v[2051] + v[1591] * v[2052] + v[1585] * v[2054]
					+ v[1600] * v[2056] + v[1604] * v[2058] + v[1602] * v[2059] + v[2060] + v[1579] * v[2063] + v[1587] * v[2064]
					+ v[1598] * v[2065] + v[2066] + v[1583] * v[2067] + v[1590] * v[2068] + v[1596] * v[2071] + v[1522] * v[2150]
					+ v[1520] * v[2151] + v[1518] * v[2153] + v[1514] * v[2157] + v[1511] * v[2158] + v[1510] * v[2161] + v[1503] * v[2165]
					+ v[1506] * v[2167] + v[1502] * v[2168] + v[1488] * v[2249] + v[1491] * v[2250] + v[1493] * v[2251] + v[1577] * v[2252]
					+ v[1581] * v[2253] + v[1582] * v[2254] + v[2044] * v[2258] + v[2045] * v[2259] + v[2049] * v[2260] + v[2050] * v[2261]
					+ v[2061] * v[2262] + v[2062] * v[2263] + v[1996] * v[2939] + v[2000] * v[2940] + v[2004] * v[2941]
					- 2e0*v[1497] * v[3177] * v[3178]) - v[2210] * v[399] - v[2211] * v[402] - v[2206] * v[404] + v[2169] * v[414]
				+ 2e0*v[2166] * v[419] + 2e0*v[2164] * v[425] + 2e0*v[2163] * v[429] + v[2160] * v[433] + 2e0*v[2156] * v[438]
				+ 2e0*v[2154] * v[442] + 2e0*v[2152] * v[446] + v[2149] * v[451] + 2e0*v[4967 + i1143])*v[962];
		v[2960] = v[2274] + (2e0*v[1367] * v[1477] - v[2149] * v[218] + v[2196] * v[411] - v[2199] * v[412] + v[2187] * v[413])
			/ 2e0;
		v[2958] = (v[2204] + v[2208]) / 2e0;
		v[2266] = v[2207] + v[2209];
		v[2267] = v[2193] + v[2205];
		v[5004] = 0e0;
		v[5005] = 0e0;
		v[5006] = 0e0;
		v[5007] = 2e0*v[2270];
		v[5008] = v[2271];
		v[5009] = v[2272];
		v[5010] = 0e0;
		v[5011] = 0e0;
		v[5012] = 0e0;
		v[5013] = 0e0;
		v[5014] = 0e0;
		v[5015] = 0e0;
		v[5016] = 0e0;
		v[5017] = 0e0;
		v[5018] = 0e0;
		v[5019] = v[2271];
		v[5020] = 2e0*v[2275];
		v[5021] = v[2276];
		v[5022] = 0e0;
		v[5023] = 0e0;
		v[5024] = 0e0;
		v[5025] = 0e0;
		v[5026] = 0e0;
		v[5027] = 0e0;
		v[5028] = 0e0;
		v[5029] = 0e0;
		v[5030] = 0e0;
		v[5031] = v[2272];
		v[5032] = v[2276];
		v[5033] = 2e0*v[2277];
		v[5034] = 0e0;
		v[5035] = 0e0;
		v[5036] = 0e0;
		v[5037] = 0e0;
		v[5038] = 0e0;
		v[5039] = 0e0;
		v[5040] = 0e0;
		v[5041] = 0e0;
		v[5042] = 0e0;
		v[5043] = 0e0;
		v[5044] = 0e0;
		v[5045] = 0e0;
		v[5046] = 0e0;
		v[5047] = 0e0;
		v[5048] = 0e0;
		v[5049] = 2e0*v[2278];
		v[5050] = v[2279];
		v[5051] = v[2280];
		v[5052] = 0e0;
		v[5053] = 0e0;
		v[5054] = 0e0;
		v[5055] = 0e0;
		v[5056] = 0e0;
		v[5057] = 0e0;
		v[5058] = 0e0;
		v[5059] = 0e0;
		v[5060] = 0e0;
		v[5061] = v[2279];
		v[5062] = 2e0*v[2283];
		v[5063] = v[2284];
		v[5064] = 0e0;
		v[5065] = 0e0;
		v[5066] = 0e0;
		v[5067] = 0e0;
		v[5068] = 0e0;
		v[5069] = 0e0;
		v[5070] = 0e0;
		v[5071] = 0e0;
		v[5072] = 0e0;
		v[5073] = v[2280];
		v[5074] = v[2284];
		v[5075] = 2e0*v[2285];
		v[5076] = v[1941];
		v[5077] = v[1949];
		v[5078] = v[1957];
		v[5079] = v[1391] * v[1477] - v[2207] + v[2209] + v[2192] * v[2957] + alphaA[2] * v[2958] + 2e0*(v[1313] * v[2273]
			+ v[2210] * v[2957] + v[2241] * v[2959]) + 2e0*alphaA[0] * v[2960] + v[2267] * v[400] + v[5003 + i1143];
		v[5080] = -(v[1392] * v[1477]) + v[2204] - v[2208] + (alphaA[2] * v[2266]) / 2e0 + (alphaA[0] * v[2267]) / 2e0
			- v[2188] * v[2957] + 2e0*alphaA[1] * (-v[2242] + v[2245] + v[2960]) + v[5015 + i1143] + 2e0*(-(v[1328] * v[2273])
				+ v[2211] * v[2957] - 4e0*v[2244] * v[962]);
		v[5081] = v[1393] * v[1477] - v[2193] + v[2205] + 2e0*alphaA[2] * (-v[2242] + v[2274]) + v[2203] * v[2957]
			+ alphaA[0] * v[2958] + 2e0*(v[1324] * v[2273] + v[2206] * v[2957] + v[2246] * v[2959]) + v[2266] * v[400] + v[5027
			+ i1143];
		v[5082] = v[1937];
		v[5083] = v[1945];
		v[5084] = v[1953];
		v[5085] = v[1353] * v[1528] - v[2138] + v[2140] + v[2123] * v[2961] + alphaB[2] * v[2962] + 2e0*(v[1286] * v[2281]
			+ v[2141] * v[2961] + v[2212] * v[2963]) + 2e0*alphaB[0] * v[2964] + v[2238] * v[406] + v[5039 + i1143];
		v[5086] = -(v[1354] * v[1528]) + v[2135] - v[2139] + (alphaB[2] * v[2237]) / 2e0 + (alphaB[0] * v[2238]) / 2e0
			- v[2119] * v[2961] + 2e0*alphaB[1] * (-v[2213] + v[2216] + v[2964]) + v[5051 + i1143] + 2e0*(-(v[1301] * v[2281])
				+ v[2142] * v[2961] - 4e0*v[2215] * v[977]);
		v[5087] = v[1355] * v[1528] - v[2124] + v[2136] + 2e0*alphaB[2] * (-v[2213] + v[2282]) + v[2134] * v[2961]
			+ alphaB[0] * v[2962] + 2e0*(v[1297] * v[2281] + v[2137] * v[2961] + v[2217] * v[2963]) + v[2237] * v[406] + v[5063
			+ i1143];
		v[2268] = ddGAp[0][0][1] * v[2143] + ddGAp[1][0][1] * v[2144] + ddGAp[2][0][1] * v[2145]
			+ ddGAp[0][1][1] * v[2146] + ddGAp[1][1][1] * v[2147] + ddGAp[2][1][1] * v[2148] + dGAp[0][1] * v[2200]
			+ dGAp[1][1] * v[2201] + dGAp[2][1] * v[2202];
		v[2269] = ddGAp[0][0][0] * v[2143] + ddGAp[1][0][0] * v[2144] + ddGAp[2][0][0] * v[2145]
			+ ddGAp[0][1][0] * v[2146] + ddGAp[1][1][0] * v[2147] + ddGAp[2][1][0] * v[2148] + dGAp[0][0] * v[2200]
			+ dGAp[1][0] * v[2201] + dGAp[2][0] * v[2202];
		Rc[i1143 - 1] += v[4207 + i1143] + (*a4)*v[4219 + i1143];
		for (i1429 = 1; i1429 <= 12; i1429++) {
			Kc[i1143 - 1][i1429 - 1] += v[2239] * v[3707 + i1429] + v[2240] * v[3719 + i1429] + v[2269] * v[3731 + i1429]
				+ v[2268] * v[3743 + i1429] + v[5075 + i1429] + (*a4)*v[5087 + i1429];
		};/* end for */
	};/* end for */
	v[2304] = v[763] * v[852];
	v[2305] = v[763] * v[854];
	v[2306] = v[763] * v[856];
	v[2307] = v[763] * v[864];
	v[2308] = v[763] * v[866];
	v[2309] = v[763] * v[868];
	v[2310] = v[762] * v[852];
	v[2311] = v[762] * v[854];
	v[2312] = v[762] * v[856];
	v[2313] = v[762] * v[864];
	v[2314] = v[762] * v[866];
	v[2315] = v[762] * v[868];
	v[2316] = v[761] * v[852];
	v[2317] = v[761] * v[854];
	v[2318] = v[761] * v[856];
	v[2319] = v[761] * v[864];
	v[2320] = v[761] * v[866];
	v[2321] = v[761] * v[868];
	v[2322] = v[237] * v[761] + v[240] * v[762] + v[243] * v[763];
	v[2323] = v[238] * v[761] + v[241] * v[762] + v[244] * v[763];
	v[2324] = v[239] * v[761] + v[242] * v[762] + v[245] * v[763];
	v[3195] = ddGAp[0][0][0] * v[2322] + ddGAp[1][0][0] * v[2323] + ddGAp[2][0][0] * v[2324];
	v[3194] = ddGAp[0][1][0] * v[2322] + ddGAp[1][1][0] * v[2323] + ddGAp[2][1][0] * v[2324];
	v[3193] = ddGAp[0][0][1] * v[2322] + ddGAp[1][0][1] * v[2323] + ddGAp[2][0][1] * v[2324];
	v[3192] = ddGAp[0][1][1] * v[2322] + ddGAp[1][1][1] * v[2323] + ddGAp[2][1][1] * v[2324];
	v[2325] = -(v[277] * v[761]) - v[280] * v[762] - v[283] * v[763];
	v[2326] = -(v[279] * v[761]) - v[282] * v[762] - v[285] * v[763];
	v[2327] = -(v[278] * v[761]) - v[281] * v[762] - v[284] * v[763];
	v[3191] = ddGBp[0][0][0] * v[2325] + ddGBp[2][0][0] * v[2326] + ddGBp[1][0][0] * v[2327];
	v[3190] = ddGBp[0][1][0] * v[2325] + ddGBp[2][1][0] * v[2326] + ddGBp[1][1][0] * v[2327];
	v[3189] = ddGBp[0][0][1] * v[2325] + ddGBp[2][0][1] * v[2326] + ddGBp[1][0][1] * v[2327];
	v[3188] = ddGBp[0][1][1] * v[2325] + ddGBp[2][1][1] * v[2326] + ddGBp[1][1][1] * v[2327];
	v[2328] = (v[2304] * v[258]) / 2e0;
	v[2329] = v[2305] * v[258];
	v[2330] = v[2306] * v[258];
	v[2331] = v[2310] * v[258];
	v[2332] = (v[2311] * v[258]) / 2e0;
	v[2333] = v[2312] * v[258];
	v[2334] = v[2316] * v[258];
	v[2335] = v[2317] * v[258];
	v[2336] = (v[2318] * v[513]) / 2e0 + v[2317] * v[518] + v[2316] * v[524] + v[2312] * v[528] + (v[2311] * v[532]) / 2e0
		+ v[2310] * v[537] + v[2306] * v[541] + v[2305] * v[545] + (v[2304] * v[550]) / 2e0;
	v[2691] = -v[2328] - v[2332] - 4e0*v[2336] * v[977];
	v[2697] = v[2328] - (v[2318] * v[258]) / 2e0 + v[2691];
	v[2695] = -v[2328] + v[2332] + v[2697];
	v[2338] = (v[218] * v[2307]) / 2e0;
	v[2339] = v[218] * v[2308];
	v[2340] = v[218] * v[2309];
	v[2341] = v[218] * v[2313];
	v[2342] = (v[218] * v[2314]) / 2e0;
	v[2343] = v[218] * v[2315];
	v[2344] = v[218] * v[2319];
	v[2345] = v[218] * v[2320];
	v[2346] = (v[2321] * v[414]) / 2e0 + v[2320] * v[419] + v[2319] * v[425] + v[2315] * v[429] + (v[2314] * v[433]) / 2e0
		+ v[2313] * v[438] + v[2309] * v[442] + v[2308] * v[446] + (v[2307] * v[451]) / 2e0;
	v[2684] = -v[2338] - v[2342] - 4e0*v[2346] * v[962];
	v[2690] = -(v[218] * v[2321]) / 2e0 + v[2338] + v[2684];
	v[2688] = -v[2338] + v[2342] + v[2690];
	v[2693] = (v[2330] + v[2334]) / 2e0;
	v[2349] = v[2329] + v[2331];
	v[2696] = v[2349] / 2e0;
	v[2350] = v[2333] + v[2335];
	v[2692] = v[2350] / 2e0;
	v[2351] = dGAp[0][1] * v[2322] + dGAp[1][1] * v[2323] + dGAp[2][1] * v[2324];
	v[2352] = dGAp[0][0] * v[2322] + dGAp[1][0] * v[2323] + dGAp[2][0] * v[2324];
	v[2353] = dGBp[0][0] * v[2325] + dGBp[2][0] * v[2326] + dGBp[1][0] * v[2327];
	v[2354] = dGBp[0][1] * v[2325] + dGBp[2][1] * v[2326] + dGBp[1][1] * v[2327];
	v[2686] = (v[2340] + v[2344]) / 2e0;
	v[2356] = v[2339] + v[2341];
	v[2689] = v[2356] / 2e0;
	v[2357] = v[2343] + v[2345];
	v[5104] = v[761];
	v[5105] = v[762];
	v[5106] = v[763];
	v[5107] = v[2339] - v[2341] + 2e0*alphaA[0] * v[2684] + alphaA[2] * v[2686] + v[2357] * v[400];
	v[5108] = -v[2340] + v[2344] + (alphaA[2] * v[2356]) / 2e0 + (alphaA[0] * v[2357]) / 2e0 + 2e0*alphaA[1] * v[2688];
	v[5109] = v[2343] - v[2345] + alphaA[0] * v[2686] + 2e0*alphaA[2] * v[2690] + v[2356] * v[400];
	v[5110] = -v[761];
	v[5111] = -v[762];
	v[5112] = -v[763];
	v[5113] = v[2329] - v[2331] + 2e0*alphaB[0] * v[2691] + alphaB[2] * v[2693] + v[2350] * v[406];
	v[5114] = -v[2330] + v[2334] + (alphaB[2] * v[2349]) / 2e0 + (alphaB[0] * v[2350]) / 2e0 + 2e0*alphaB[1] * v[2695];
	v[5115] = v[2333] - v[2335] + alphaB[0] * v[2693] + 2e0*alphaB[2] * v[2697] + v[2349] * v[406];
	v[2685] = v[2357] / 2e0;
	for (i2302 = 1; i2302 <= 12; i2302++) {
		i2997 = (i2302 == 11 ? 1 : 0);
		i2996 = (i2302 == 10 ? 1 : 0);
		i2995 = (i2302 == 12 ? 1 : 0);
		i2994 = (i2302 == 5 ? 1 : 0);
		i2993 = (i2302 == 4 ? 1 : 0);
		i2992 = (i2302 == 6 ? 1 : 0);
		v[2388] = v[3731 + i2302];
		v[2387] = v[3743 + i2302];
		v[2383] = v[3719 + i2302];
		v[2382] = v[3707 + i2302];
		v[2364] = v[3783 + i2302];
		v[2365] = v[3807 + i2302];
		v[2366] = v[3795 + i2302];
		v[2367] = v[3831 + i2302];
		v[2368] = v[3855 + i2302];
		v[2369] = v[3843 + i2302];
		v[2371] = v[4235 + i2302];
		v[2372] = v[3771 + i2302];
		v[2413] = -4e0*v[2372] * v[962];
		v[2373] = v[5119 + i2302];
		v[2375] = v[4259 + i2302];
		v[2377] = v[4307 + i2302];
		v[2378] = v[3819 + i2302];
		v[2424] = -4e0*v[2378] * v[977];
		v[2379] = v[5131 + i2302];
		v[2381] = v[4331 + i2302];
		v[2384] = dGBp[1][1] * v[2382] + dGBp[1][0] * v[2383];
		v[2385] = dGBp[2][1] * v[2382] + dGBp[2][0] * v[2383];
		v[2386] = dGBp[0][1] * v[2382] + dGBp[0][0] * v[2383];
		v[2389] = dGAp[2][1] * v[2387] + dGAp[2][0] * v[2388];
		v[2390] = dGAp[1][1] * v[2387] + dGAp[1][0] * v[2388];
		v[2391] = dGAp[0][1] * v[2387] + dGAp[0][0] * v[2388];
		v[2392] = -i2992 + v[2364];
		v[2394] = i2992 + v[2364];
		v[2395] = -i2993 + v[2365];
		v[2397] = i2993 + v[2365];
		v[2398] = i2994 + v[2366];
		v[2400] = -i2994 + v[2366];
		v[2401] = -i2995 + v[2367];
		v[2403] = i2995 + v[2367];
		v[2404] = -i2996 + v[2368];
		v[2406] = i2996 + v[2368];
		v[2407] = i2997 + v[2369];
		v[2409] = -i2997 + v[2369];
		v[2411] = (v[218] * v[2371]) / 2e0 + v[2372] * v[2410];
		v[2412] = v[218] * v[2392] + v[2413] * v[419];
		v[2414] = v[218] * v[2398] + v[2413] * v[425];
		v[2415] = v[218] * v[2394] + v[2413] * v[429];
		v[2416] = (v[218] * v[2373] + v[2413] * v[433]) / 2e0;
		v[2417] = v[218] * v[2395] + v[2413] * v[438];
		v[2418] = v[218] * v[2400] + v[2413] * v[442];
		v[2419] = v[218] * v[2397] + v[2413] * v[446];
		v[2420] = (v[218] * v[2375] + v[2413] * v[451]) / 2e0;
		v[2422] = v[2378] * v[2421] + (v[2377] * v[258]) / 2e0;
		v[2423] = v[2401] * v[258] + v[2424] * v[518];
		v[2425] = v[2407] * v[258] + v[2424] * v[524];
		v[2426] = v[2403] * v[258] + v[2424] * v[528];
		v[2427] = (v[2379] * v[258] + v[2424] * v[532]) / 2e0;
		v[2428] = v[2404] * v[258] + v[2424] * v[537];
		v[2429] = v[2409] * v[258] + v[2424] * v[541];
		v[2430] = v[2406] * v[258] + v[2424] * v[545];
		v[2431] = (v[2381] * v[258] + v[2424] * v[550]) / 2e0;
		v[2433] = v[2389] * v[239] + v[238] * v[2390] + v[237] * v[2391] - v[2386] * v[277] - v[2384] * v[278] - v[2385] * v[279]
			+ v[4107 + i2302] + v[2425] * v[852] + v[2423] * v[854] + v[2422] * v[856] + v[2414] * v[864] + v[2412] * v[866]
			+ v[2411] * v[868];
		v[2435] = v[2391] * v[240] + v[2390] * v[241] + v[2389] * v[242] - v[2386] * v[280] - v[2384] * v[281] - v[2385] * v[282]
			+ v[4095 + i2302] + v[2428] * v[852] + v[2427] * v[854] + v[2426] * v[856] + v[2417] * v[864] + v[2416] * v[866]
			+ v[2415] * v[868];
		v[2436] = v[2411] * v[761] + v[2415] * v[762] + v[2418] * v[763];
		v[2437] = v[2412] * v[761] + v[2416] * v[762] + v[2419] * v[763];
		v[2438] = v[2414] * v[761] + v[2417] * v[762] + v[2420] * v[763];
		v[2439] = v[2422] * v[761] + v[2426] * v[762] + v[2429] * v[763];
		v[2440] = v[2423] * v[761] + v[2427] * v[762] + v[2430] * v[763];
		v[2442] = v[2391] * v[243] + v[2390] * v[244] + v[2389] * v[245] - v[2386] * v[283] - v[2384] * v[284] - v[2385] * v[285]
			+ v[4083 + i2302] + v[2431] * v[852] + v[2430] * v[854] + v[2429] * v[856] + v[2420] * v[864] + v[2419] * v[866]
			+ v[2418] * v[868];
		v[2699] = v[2288] * v[2433] - v[1909] * v[2442] + (*cn)*v[2435] * v[751];
		v[2698] = v[2288] * v[2435] + v[2290] * v[2442] + (*cn)*v[2433] * v[749];
		v[2700] = v[2290] * v[2433] - v[1909] * v[2435] + (*cn)*v[2442] * v[753];
		v[2443] = v[2425] * v[761] + v[2428] * v[762] + v[2431] * v[763];
		v[2444] = (*cn)*(v[1981] * v[2433] + v[1982] * v[2435] + v[1983] * v[2442]);
		v[2553] = v[1426] * v[2444];
		v[2539] = v[2298] * v[2444];
		v[2445] = (*cn)*(v[1985] * v[2433] + v[1986] * v[2435] + v[1987] * v[2442]);
		v[2546] = v[2297] * v[2445];
		v[3006] = v[1425] * v[2444] + v[2546];
		v[3005] = v[1423] * v[2445] + v[2539];
		v[2446] = (*cn)*(v[1989] * v[2433] + v[1990] * v[2435] + v[1991] * v[2442]);
		v[2551] = v[2296] * v[2446];
		v[3007] = v[1424] * v[2445] + v[2551];
		v[2548] = v[2446] * v[2715];
		v[2542] = v[2446] * v[2717];
		v[2447] = (*cn)*(v[1993] * v[2433] + v[1994] * v[2435] + v[1995] * v[2442]);
		v[2611] = v[1420] * v[2447];
		v[2597] = v[2295] * v[2447];
		v[2448] = (*cn)*(v[1997] * v[2433] + v[1998] * v[2435] + v[1999] * v[2442]);
		v[2604] = v[2294] * v[2448];
		v[3009] = v[1419] * v[2447] + v[2604];
		v[3008] = v[1417] * v[2448] + v[2597];
		v[2449] = (*cn)*(v[2001] * v[2433] + v[2002] * v[2435] + v[2003] * v[2442]);
		v[2609] = v[2293] * v[2449];
		v[5228] = v[2698];
		v[5229] = v[2699];
		v[5230] = v[2700];
		v[5231] = v[2609] + v[2448] * v[2710] + v[2447] * v[2712];
		v[5232] = v[1417] * v[2447] + v[1418] * v[2449] + v[2604];
		v[5233] = v[1419] * v[2448] + v[1420] * v[2449] + v[2597];
		v[5234] = -v[2698];
		v[5235] = -v[2699];
		v[5236] = -v[2700];
		v[5237] = v[2551] + v[2445] * v[2715] + v[2444] * v[2717];
		v[5238] = v[1423] * v[2444] + v[1424] * v[2446] + v[2546];
		v[5239] = v[1425] * v[2445] + v[1426] * v[2446] + v[2539];
		v[3010] = v[1418] * v[2448] + v[2609];
		v[2606] = v[2449] * v[2710];
		v[2600] = v[2449] * v[2712];
		v[2450] = -(v[2698] * v[345]);
		v[2452] = -(v[2698] * v[355]);
		v[2453] = -(v[2698] * v[335]);
		v[2454] = -(v[2699] * v[335]);
		v[2456] = -(v[2699] * v[355]);
		v[2457] = -(v[2700] * v[355]);
		v[2458] = -(v[2699] * v[345]);
		v[2460] = -(v[2700] * v[335]);
		v[2461] = -(v[2700] * v[345]);
		v[2462] = v[2698] * v[319];
		v[2463] = v[2698] * v[329];
		v[2464] = v[2698] * v[309];
		v[2465] = v[2699] * v[309];
		v[2466] = v[2699] * v[329];
		v[2467] = v[2700] * v[329];
		v[2468] = v[2699] * v[319];
		v[2469] = v[2700] * v[309];
		v[2470] = v[2700] * v[319];
		v[2475] = (*cn)*(v[2433] * v[2471] + v[2435] * v[2472] + v[2442] * v[3185]);
		v[2480] = (*cn)*(v[2433] * v[2476] + v[2442] * v[2479] + v[2435] * v[3186]);
		v[2485] = (*cn)*(v[2435] * v[2483] + v[2442] * v[2484] + v[2433] * v[3187]);
		v[2998] = ((v[2485] * v[356] + v[2480] * v[357] + v[2475] * v[358])*v[369]) / v[1264];
		v[2487] = v[2998] * v[358] + v[2475] * v[370];
		v[2488] = v[2998] * v[357] + v[2480] * v[370];
		v[2489] = v[2998] * v[356] + v[2485] * v[370];
		v[2490] = -(GBp[2] * v[2487]) - v[2385] * v[763];
		v[2491] = -(GBp[1] * v[2487]) - v[2384] * v[763];
		v[2492] = -(GBp[0] * v[2487]) - v[2386] * v[763];
		v[2493] = GAp[2] * v[2487] + v[2389] * v[763];
		v[2494] = GAp[1] * v[2487] + v[2390] * v[763];
		v[2495] = GAp[0] * v[2487] + v[2391] * v[763];
		v[2496] = -(GBp[2] * v[2488]) - v[2385] * v[762];
		v[2497] = -(GBp[1] * v[2488]) - v[2384] * v[762];
		v[2498] = -(GBp[0] * v[2488]) - v[2386] * v[762];
		v[2499] = GAp[2] * v[2488] + v[2389] * v[762];
		v[2500] = GAp[1] * v[2488] + v[2390] * v[762];
		v[2501] = GAp[0] * v[2488] + v[2391] * v[762];
		v[2502] = -(GBp[2] * v[2489]) - v[2385] * v[761];
		v[2503] = -(GBp[1] * v[2489]) - v[2384] * v[761];
		v[2504] = -(GBp[0] * v[2489]) - v[2386] * v[761];
		v[2505] = GAp[2] * v[2489] + v[2389] * v[761];
		v[2506] = GAp[1] * v[2489] + v[2390] * v[761];
		v[2507] = GAp[0] * v[2489] + v[2391] * v[761];
		v[2508] = v[1785] * v[2444];
		v[2509] = v[1784] * v[2444];
		v[2510] = v[1783] * v[2444];
		v[2511] = v[1780] * v[2445];
		v[2512] = v[2444] * v[342] + v[2445] * v[344];
		v[2513] = v[1778] * v[2445];
		v[2514] = v[2508] + v[2511];
		v[2515] = v[1779] * v[2445] + v[2512] * v[2999];
		v[2516] = v[1773] * v[2446];
		v[2517] = v[2444] * v[336] + v[2446] * v[344];
		v[2518] = v[2445] * v[336] + v[2446] * v[342];
		v[2519] = v[1775] * v[2446] + v[2517] * v[3000];
		v[2520] = v[2515] + v[2519];
		v[2521] = v[1774] * v[2446] + v[2518] * v[3001];
		v[2522] = v[2509] + v[2521];
		v[2523] = v[1770] * v[2447];
		v[2524] = v[1769] * v[2447];
		v[2525] = v[1768] * v[2447];
		v[2526] = v[1765] * v[2448];
		v[2527] = v[2447] * v[316] + v[2448] * v[318];
		v[2528] = v[1763] * v[2448];
		v[2529] = v[2523] + v[2526];
		v[2530] = v[1764] * v[2448] + v[2527] * v[3002];
		v[2531] = v[1758] * v[2449];
		v[2532] = v[2447] * v[310] + v[2449] * v[318];
		v[2533] = v[2448] * v[310] + v[2449] * v[316];
		v[2534] = v[1760] * v[2449] + v[2532] * v[3003];
		v[2535] = v[2530] + v[2534];
		v[2536] = v[1759] * v[2449] + v[2533] * v[3004];
		v[2537] = v[2524] + v[2536];
		v[2538] = QBi[2][2] * v[2490] + QBi[2][1] * v[2491] + QBi[2][0] * v[2492] + (v[2542] + v[3005])*v[344];
		v[2541] = QBi[1][2] * v[2490] + QBi[1][1] * v[2491] + QBi[1][0] * v[2492] + v[2518] * v[2717] + v[3005] * v[342];
		v[2543] = QBi[0][2] * v[2490] + QBi[0][1] * v[2491] + QBi[0][0] * v[2492] + (v[2539] + v[2542])*v[336];
		v[2544] = QBi[2][2] * v[2496] + QBi[2][1] * v[2497] + QBi[2][0] * v[2498] + v[2517] * v[2715] + v[3006] * v[344];
		v[2547] = QBi[1][2] * v[2496] + QBi[1][1] * v[2497] + QBi[1][0] * v[2498] + (v[2548] + v[3006])*v[342];
		v[2549] = QBi[0][2] * v[2496] + QBi[0][1] * v[2497] + QBi[0][0] * v[2498] + (v[2546] + v[2548])*v[336];
		v[2550] = QBi[2][2] * v[2502] + QBi[2][1] * v[2503] + QBi[2][0] * v[2504] + v[1424] * v[2512] + (v[2551] + v[2553]
			)*v[344];
		v[2552] = QBi[1][2] * v[2502] + QBi[1][1] * v[2503] + QBi[1][0] * v[2504] + v[3007] * v[342];
		v[2555] = QBi[0][2] * v[2502] + QBi[0][1] * v[2503] + QBi[0][0] * v[2504] + (v[2553] + v[3007])*v[336];
		v[2556] = -(v[2450] * v[856]);
		v[2557] = -(v[2450] * v[854]);
		v[2558] = -(v[2450] * v[852]);
		v[2559] = -(v[2452] * v[856]);
		v[2560] = -(v[2452] * v[852]);
		v[2561] = -(v[2452] * v[854]);
		v[2562] = -(v[2453] * v[854]);
		v[2563] = -(v[2453] * v[852]);
		v[2564] = -(v[2453] * v[856]);
		v[2565] = -(v[2454] * v[856]);
		v[2566] = -(v[2454] * v[854]);
		v[2567] = -(v[2454] * v[852]);
		v[2568] = -(v[2456] * v[852]);
		v[2569] = -(v[2456] * v[856]);
		v[2570] = -(v[2456] * v[854]);
		v[2571] = -(v[2457] * v[854]);
		v[2572] = -(v[2457] * v[856]);
		v[2573] = -(v[2457] * v[852]);
		v[2574] = v[2562] + v[2565] + v[2568] + v[2571] + 2e0*v[2513] * v[517] - v[2520] * v[520] - v[2514] * v[523];
		v[2575] = -(v[2458] * v[856]);
		v[2576] = -(v[2458] * v[852]);
		v[2577] = -(v[2458] * v[854]);
		v[2578] = -v[2557] - v[2560] - v[2572] - v[2575] - v[2520] * v[517] + 2e0*v[2516] * v[520] + v[2522] * v[523];
		v[2579] = v[2317] * v[2424] + v[2552] * v[258] - v[2562] * v[510] + v[2557] * v[511] - v[2561] * v[512];
		v[2580] = -(v[2460] * v[856]);
		v[2581] = -(v[2460] * v[854]);
		v[2582] = -(v[2460] * v[852]);
		v[2583] = -(v[2461] * v[854]);
		v[2584] = -(v[2461] * v[856]);
		v[2585] = -(v[2461] * v[852]);
		v[2586] = -v[2563] - v[2576] - v[2580] - v[2583] - v[2514] * v[517] + v[2522] * v[520] + 2e0*v[2510] * v[523];
		v[2587] = v[2316] * v[2424] + v[2550] * v[258] - v[2563] * v[510] + v[2558] * v[511] - v[2560] * v[512];
		v[2588] = v[2312] * v[2424] + v[2549] * v[258] - v[2565] * v[510] + v[2575] * v[511] - v[2569] * v[512];
		v[2589] = v[2559] + v[2570];
		v[2590] = v[2310] * v[2424] + v[2544] * v[258] - v[2567] * v[510] + v[2576] * v[511] - v[2568] * v[512];
		v[2591] = v[2306] * v[2424] + v[2543] * v[258] - v[2580] * v[510] + v[2584] * v[511] - v[2572] * v[512];
		v[2592] = v[2305] * v[2424] + v[2541] * v[258] - v[2581] * v[510] + v[2583] * v[511] - v[2571] * v[512];
		v[2593] = v[2566] + v[2582];
		v[2594] = v[2556] + v[2585];
		v[2596] = QAi[2][2] * v[2493] + QAi[2][1] * v[2494] + QAi[2][0] * v[2495] + (v[2600] + v[3008])*v[318];
		v[2599] = QAi[1][2] * v[2493] + QAi[1][1] * v[2494] + QAi[1][0] * v[2495] + v[2533] * v[2712] + v[3008] * v[316];
		v[2601] = QAi[0][2] * v[2493] + QAi[0][1] * v[2494] + QAi[0][0] * v[2495] + (v[2597] + v[2600])*v[310];
		v[2602] = QAi[2][2] * v[2499] + QAi[2][1] * v[2500] + QAi[2][0] * v[2501] + v[2532] * v[2710] + v[3009] * v[318];
		v[2605] = QAi[1][2] * v[2499] + QAi[1][1] * v[2500] + QAi[1][0] * v[2501] + (v[2606] + v[3009])*v[316];
		v[2607] = QAi[0][2] * v[2499] + QAi[0][1] * v[2500] + QAi[0][0] * v[2501] + (v[2604] + v[2606])*v[310];
		v[2608] = QAi[2][2] * v[2505] + QAi[2][1] * v[2506] + QAi[2][0] * v[2507] + v[1418] * v[2527] + (v[2609] + v[2611]
			)*v[318];
		v[2610] = QAi[1][2] * v[2505] + QAi[1][1] * v[2506] + QAi[1][0] * v[2507] + v[3010] * v[316];
		v[2613] = QAi[0][2] * v[2505] + QAi[0][1] * v[2506] + QAi[0][0] * v[2507] + (v[2611] + v[3010])*v[310];
		v[2614] = v[2462] * v[868];
		v[2615] = v[2462] * v[866];
		v[2616] = v[2462] * v[864];
		v[2617] = v[2463] * v[868];
		v[2618] = v[2463] * v[864];
		v[2619] = v[2463] * v[866];
		v[2620] = v[2464] * v[866];
		v[2621] = v[2464] * v[864];
		v[2622] = v[2464] * v[868];
		v[2623] = v[2465] * v[868];
		v[2624] = v[2465] * v[866];
		v[2625] = v[2465] * v[864];
		v[2626] = v[2466] * v[864];
		v[2627] = v[2466] * v[868];
		v[2628] = v[2466] * v[866];
		v[2629] = v[2467] * v[866];
		v[2630] = v[2467] * v[868];
		v[2631] = v[2467] * v[864];
		v[2632] = v[2620] + v[2623] + v[2626] + v[2629] + 2e0*v[2528] * v[418] - v[2535] * v[421] - v[2529] * v[424];
		v[2633] = v[2468] * v[868];
		v[2634] = v[2468] * v[864];
		v[2635] = v[2468] * v[866];
		v[2636] = -v[2615] - v[2618] - v[2630] - v[2633] - v[2535] * v[418] + 2e0*v[2531] * v[421] + v[2537] * v[424];
		v[2637] = v[2320] * v[2413] + v[218] * v[2610] - v[2620] * v[411] + v[2615] * v[412] - v[2619] * v[413];
		v[2638] = v[2469] * v[868];
		v[2639] = v[2469] * v[866];
		v[2640] = v[2469] * v[864];
		v[2641] = v[2470] * v[866];
		v[2642] = v[2470] * v[868];
		v[2643] = v[2470] * v[864];
		v[2644] = -v[2621] - v[2634] - v[2638] - v[2641] - v[2529] * v[418] + v[2537] * v[421] + 2e0*v[2525] * v[424];
		v[2645] = v[2319] * v[2413] + v[218] * v[2608] - v[2621] * v[411] + v[2616] * v[412] - v[2618] * v[413];
		v[2646] = v[2315] * v[2413] + v[218] * v[2607] - v[2623] * v[411] + v[2633] * v[412] - v[2627] * v[413];
		v[2647] = v[2617] + v[2628];
		v[2648] = v[2313] * v[2413] + v[218] * v[2602] - v[2625] * v[411] + v[2634] * v[412] - v[2626] * v[413];
		v[2649] = v[2309] * v[2413] + v[218] * v[2601] - v[2638] * v[411] + v[2642] * v[412] - v[2630] * v[413];
		v[2650] = v[2308] * v[2413] + v[218] * v[2599] - v[2639] * v[411] + v[2641] * v[412] - v[2629] * v[413];
		v[2651] = v[2624] + v[2640];
		v[2652] = v[2614] + v[2643];
		v[2654] = (-2e0*v[2508] + 2e0*v[2511] - v[2564] * v[513] - 2e0*v[2562] * v[518] - 2e0*v[2563] * v[524]
			- 2e0*v[2565] * v[528] - v[2566] * v[532] - 2e0*v[2567] * v[537] - 2e0*v[2580] * v[541] - 2e0*v[2581] * v[545]
			- v[2582] * v[550]) / 2e0;
		v[2655] = (v[2318] * v[2424] + v[2555] * v[258] - v[2564] * v[510] + v[2556] * v[511] - v[2559] * v[512]) / 2e0;
		v[2657] = -v[2509] + v[2521] + (v[2556] * v[513]) / 2e0 + v[2557] * v[518] + v[2558] * v[524] + v[2575] * v[528] +
			(v[2577] * v[532]) / 2e0 + v[2576] * v[537] + v[2584] * v[541] + v[2583] * v[545] + (v[2585] * v[550]) / 2e0;
		v[2658] = (v[2311] * v[2424] + v[2547] * v[258] - v[2566] * v[510] + v[2577] * v[511] - v[2570] * v[512]) / 2e0;
		v[2659] = (-2e0*v[2515] + 2e0*v[2519] - v[2559] * v[513] - 2e0*v[2561] * v[518] - 2e0*v[2560] * v[524]
			- 2e0*v[2569] * v[528] - v[2570] * v[532] - 2e0*v[2568] * v[537] - 2e0*v[2572] * v[541] - 2e0*v[2571] * v[545]
			- v[2573] * v[550]) / 2e0;
		v[2694] = v[1531] * v[2657] - v[2658] + v[1533] * v[2659] + 8e0*v[1413] * (v[2336] * v[2378] - v[2654] * v[405]) - 2e0*
			(v[2318] * v[2377] + v[2311] * v[2379] + v[2304] * v[2381] + 2e0*v[2317] * v[2401] + 2e0*v[2312] * v[2403]
				+ 2e0*v[2310] * v[2404] + 2e0*v[2305] * v[2406] + 2e0*v[2316] * v[2407] + 2e0*v[2306] * v[2409] + 2e0*v[2558]
				- 2e0*v[2561] - 2e0*v[2567] + 2e0*v[2569] + alphaB[1] * v[2574] - alphaB[0] * v[2578] + 4e0*(v[1782] * v[2444]
					+ v[1777] * v[2445] + v[1772] * v[2446] + v[2510] + v[1613] * v[2512] + v[2513] + v[2516] + v[1611] * v[2517]
					+ v[1607] * v[2518])*v[258] + 2e0*v[2581] - 2e0*v[2584] - alphaB[2] * v[2586] - v[2593] * v[405] - v[2594] * v[408]
				- v[2589] * v[410] + v[2555] * v[513] + 2e0*v[2552] * v[518] + 2e0*v[2550] * v[524] + 2e0*v[2549] * v[528]
				+ v[2547] * v[532] + 2e0*v[2544] * v[537] + 2e0*v[2543] * v[541] + 2e0*v[2541] * v[545] + v[2538] * v[550])*v[977];
		v[3014] = v[2694] + (-(v[2304] * v[2424]) - v[2538] * v[258] + v[2582] * v[510] - v[2585] * v[511] + v[2573] * v[512])
			/ 2e0;
		v[3013] = (v[2587] + v[2591]) / 2e0;
		v[2662] = v[2590] + v[2592];
		v[2663] = v[2579] + v[2588];
		v[2664] = -(QBi[0][2] * v[2439]) - QBi[1][2] * v[2440] - QBi[2][2] * v[2443] - v[2489] * v[279] - v[2488] * v[282]
			- v[2487] * v[285] + v[2453] * v[561] + v[2450] * v[562] + v[2452] * v[563] + v[2454] * v[576] + v[2458] * v[577]
			+ v[2456] * v[578] + v[2460] * v[591] + v[2461] * v[592] + v[2457] * v[593];
		v[2665] = -(QBi[0][1] * v[2439]) - QBi[1][1] * v[2440] - QBi[2][1] * v[2443] - v[2489] * v[278] - v[2488] * v[281]
			- v[2487] * v[284] + v[2453] * v[558] + v[2450] * v[559] + v[2452] * v[560] + v[2454] * v[573] + v[2458] * v[574]
			+ v[2456] * v[575] + v[2460] * v[588] + v[2461] * v[589] + v[2457] * v[590];
		v[2666] = -(QBi[0][0] * v[2439]) - QBi[1][0] * v[2440] - QBi[2][0] * v[2443] - v[2489] * v[277] - v[2488] * v[280]
			- v[2487] * v[283] + v[2453] * v[555] + v[2450] * v[556] + v[2452] * v[557] + v[2454] * v[570] + v[2458] * v[571]
			+ v[2456] * v[572] + v[2460] * v[585] + v[2461] * v[586] + v[2457] * v[587];
		v[2667] = dGBp[2][1] * v[2664] + dGBp[1][1] * v[2665] + dGBp[0][1] * v[2666] + v[2382] * v[3188] + v[2383] * v[3189];
		v[2668] = dGBp[2][0] * v[2664] + dGBp[1][0] * v[2665] + dGBp[0][0] * v[2666] + v[2382] * v[3190] + v[2383] * v[3191];
		v[2669] = (-2e0*v[2523] + 2e0*v[2526] - v[2622] * v[414] - 2e0*v[2620] * v[419] - 2e0*v[2621] * v[425]
			- 2e0*v[2623] * v[429] - v[2624] * v[433] - 2e0*v[2625] * v[438] - 2e0*v[2638] * v[442] - 2e0*v[2639] * v[446]
			- v[2640] * v[451]) / 2e0;
		v[2670] = (v[2321] * v[2413] + v[218] * v[2613] - v[2622] * v[411] + v[2614] * v[412] - v[2617] * v[413]) / 2e0;
		v[2672] = -v[2524] + v[2536] + (v[2614] * v[414]) / 2e0 + v[2615] * v[419] + v[2616] * v[425] + v[2633] * v[429] +
			(v[2635] * v[433]) / 2e0 + v[2634] * v[438] + v[2642] * v[442] + v[2641] * v[446] + (v[2643] * v[451]) / 2e0;
		v[2673] = (v[2314] * v[2413] + v[218] * v[2605] - v[2624] * v[411] + v[2635] * v[412] - v[2628] * v[413]) / 2e0;
		v[2674] = (-2e0*v[2530] + 2e0*v[2534] - v[2617] * v[414] - 2e0*v[2619] * v[419] - 2e0*v[2618] * v[425]
			- 2e0*v[2627] * v[429] - v[2628] * v[433] - 2e0*v[2626] * v[438] - 2e0*v[2630] * v[442] - 2e0*v[2629] * v[446]
			- v[2631] * v[451]) / 2e0;
		v[2687] = v[1480] * v[2672] - v[2673] + v[1482] * v[2674] + 8e0*v[1411] * (v[2346] * v[2372] - v[2669] * v[399]) - 2e0*
			(v[2321] * v[2371] + v[2314] * v[2373] + v[2307] * v[2375] + 2e0*v[2320] * v[2392] + 2e0*v[2315] * v[2394]
				+ 2e0*v[2313] * v[2395] + 2e0*v[2308] * v[2397] + 2e0*v[2319] * v[2398] + 2e0*v[2309] * v[2400] + 4e0*v[218] *
				(v[1767] * v[2447] + v[1762] * v[2448] + v[1757] * v[2449] + v[2525] + v[1585] * v[2527] + v[2528] + v[2531]
					+ v[1583] * v[2532] + v[1579] * v[2533]) + 2e0*v[2616] - 2e0*v[2619] - 2e0*v[2625] + 2e0*v[2627]
				+ alphaA[1] * v[2632] - alphaA[0] * v[2636] + 2e0*v[2639] - 2e0*v[2642] - alphaA[2] * v[2644] - v[2651] * v[399]
				- v[2652] * v[402] - v[2647] * v[404] + v[2613] * v[414] + 2e0*v[2610] * v[419] + 2e0*v[2608] * v[425]
				+ 2e0*v[2607] * v[429] + v[2605] * v[433] + 2e0*v[2602] * v[438] + 2e0*v[2601] * v[442] + 2e0*v[2599] * v[446]
				+ v[2596] * v[451])*v[962];
		v[3012] = v[2687] + (-(v[2307] * v[2413]) - v[218] * v[2596] + v[2640] * v[411] - v[2643] * v[412] + v[2631] * v[413])
			/ 2e0;
		v[3011] = (v[2645] + v[2649]) / 2e0;
		v[2677] = v[2648] + v[2650];
		v[2678] = v[2637] + v[2646];
		v[5144] = 0e0;
		v[5145] = 0e0;
		v[5146] = 0e0;
		v[5147] = 2e0*v[2684];
		v[5148] = v[2685];
		v[5149] = v[2686];
		v[5150] = 0e0;
		v[5151] = 0e0;
		v[5152] = 0e0;
		v[5153] = 0e0;
		v[5154] = 0e0;
		v[5155] = 0e0;
		v[5156] = 0e0;
		v[5157] = 0e0;
		v[5158] = 0e0;
		v[5159] = v[2685];
		v[5160] = 2e0*v[2688];
		v[5161] = v[2689];
		v[5162] = 0e0;
		v[5163] = 0e0;
		v[5164] = 0e0;
		v[5165] = 0e0;
		v[5166] = 0e0;
		v[5167] = 0e0;
		v[5168] = 0e0;
		v[5169] = 0e0;
		v[5170] = 0e0;
		v[5171] = v[2686];
		v[5172] = v[2689];
		v[5173] = 2e0*v[2690];
		v[5174] = 0e0;
		v[5175] = 0e0;
		v[5176] = 0e0;
		v[5177] = 0e0;
		v[5178] = 0e0;
		v[5179] = 0e0;
		v[5180] = 0e0;
		v[5181] = 0e0;
		v[5182] = 0e0;
		v[5183] = 0e0;
		v[5184] = 0e0;
		v[5185] = 0e0;
		v[5186] = 0e0;
		v[5187] = 0e0;
		v[5188] = 0e0;
		v[5189] = 2e0*v[2691];
		v[5190] = v[2692];
		v[5191] = v[2693];
		v[5192] = 0e0;
		v[5193] = 0e0;
		v[5194] = 0e0;
		v[5195] = 0e0;
		v[5196] = 0e0;
		v[5197] = 0e0;
		v[5198] = 0e0;
		v[5199] = 0e0;
		v[5200] = 0e0;
		v[5201] = v[2692];
		v[5202] = 2e0*v[2695];
		v[5203] = v[2696];
		v[5204] = 0e0;
		v[5205] = 0e0;
		v[5206] = 0e0;
		v[5207] = 0e0;
		v[5208] = 0e0;
		v[5209] = 0e0;
		v[5210] = 0e0;
		v[5211] = 0e0;
		v[5212] = 0e0;
		v[5213] = v[2693];
		v[5214] = v[2696];
		v[5215] = 2e0*v[2697];
		v[5216] = v[2489];
		v[5217] = v[2488];
		v[5218] = v[2487];
		v[5219] = -v[2648] + v[2650] + v[2636] * v[2957] + 2e0*(v[2651] * v[2957] + v[2669] * v[2959]) + alphaA[2] * v[3011]
			+ 2e0*alphaA[0] * v[3012] + v[2678] * v[400] + v[5143 + i2302];
		v[5220] = v[2645] - v[2649] + (alphaA[2] * v[2677]) / 2e0 + (alphaA[0] * v[2678]) / 2e0 - v[2632] * v[2957]
			+ 2e0*alphaA[1] * (-v[2670] + v[2673] + v[3012]) + v[5155 + i2302] + 2e0*(v[2652] * v[2957] - 4e0*v[2672] * v[962]);
		v[5221] = -v[2637] + v[2646] + 2e0*alphaA[2] * (-v[2670] + v[2687]) + v[2644] * v[2957] + 2e0*(v[2647] * v[2957]
			+ v[2674] * v[2959]) + alphaA[0] * v[3011] + v[2677] * v[400] + v[5167 + i2302];
		v[5222] = -v[2489];
		v[5223] = -v[2488];
		v[5224] = -v[2487];
		v[5225] = -v[2590] + v[2592] + v[2578] * v[2961] + 2e0*(v[2593] * v[2961] + v[2654] * v[2963]) + alphaB[2] * v[3013]
			+ 2e0*alphaB[0] * v[3014] + v[2663] * v[406] + v[5179 + i2302];
		v[5226] = v[2587] - v[2591] + (alphaB[2] * v[2662]) / 2e0 + (alphaB[0] * v[2663]) / 2e0 - v[2574] * v[2961]
			+ 2e0*alphaB[1] * (-v[2655] + v[2658] + v[3014]) + v[5191 + i2302] + 2e0*(v[2594] * v[2961] - 4e0*v[2657] * v[977]);
		v[5227] = -v[2579] + v[2588] + 2e0*alphaB[2] * (-v[2655] + v[2694]) + v[2586] * v[2961] + 2e0*(v[2589] * v[2961]
			+ v[2659] * v[2963]) + alphaB[0] * v[3013] + v[2662] * v[406] + v[5203 + i2302];
		v[2679] = QAi[0][0] * v[2436] + QAi[1][0] * v[2437] + QAi[2][0] * v[2438] + v[243] * v[2487] + v[240] * v[2488]
			+ v[237] * v[2489] + v[2464] * v[456] + v[2462] * v[457] + v[2463] * v[458] + v[2465] * v[471] + v[2468] * v[472]
			+ v[2466] * v[473] + v[2469] * v[486] + v[2470] * v[487] + v[2467] * v[488];
		v[2680] = QAi[0][1] * v[2436] + QAi[1][1] * v[2437] + QAi[2][1] * v[2438] + v[244] * v[2487] + v[241] * v[2488]
			+ v[238] * v[2489] + v[2464] * v[459] + v[2462] * v[460] + v[2463] * v[461] + v[2465] * v[474] + v[2468] * v[475]
			+ v[2466] * v[476] + v[2469] * v[489] + v[2470] * v[490] + v[2467] * v[491];
		v[2681] = QAi[0][2] * v[2436] + QAi[1][2] * v[2437] + QAi[2][2] * v[2438] + v[245] * v[2487] + v[242] * v[2488]
			+ v[239] * v[2489] + v[2464] * v[462] + v[2462] * v[463] + v[2463] * v[464] + v[2465] * v[477] + v[2468] * v[478]
			+ v[2466] * v[479] + v[2469] * v[492] + v[2470] * v[493] + v[2467] * v[494];
		v[2682] = dGAp[0][1] * v[2679] + dGAp[1][1] * v[2680] + dGAp[2][1] * v[2681] + v[2387] * v[3192] + v[2388] * v[3193];
		v[2683] = dGAp[0][0] * v[2679] + dGAp[1][0] * v[2680] + dGAp[2][0] * v[2681] + v[2387] * v[3194] + v[2388] * v[3195];
		Rc[i2302 - 1] += v[2354] * v[2382] + v[2353] * v[2383] + v[2351] * v[2387] + v[2352] * v[2388] + v[5103 + i2302];
		for (i2362 = 1; i2362 <= 12; i2362++) {
			Kc[i2302 - 1][i2362 - 1] += v[2667] * v[3707 + i2362] + v[2668] * v[3719 + i2362] + v[2683] * v[3731 + i2362]
				+ v[2682] * v[3743 + i2362] + v[5215 + i2362] + (*a4)*v[5227 + i2362];
		};/* end for */
	};/* end for */
#pragma endregion


}

//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
double RigidNURBS_1_RigidNURBS_1::ObjectivePhase1(Matrix& mc)
{
	//AceGen variables or pointers
	double v[2000];
	EvaluateNURBSDerivatives_p(mc);
	EvaluateNURBSDOFsVariables();
	
	double Ob;

	v[172] = Power(alphaA[0], 2);
	v[170] = 0.5e0*alphaA[0] * alphaA[1];
	v[165] = Power(alphaA[1], 2);
	v[177] = 0.5e0*alphaA[1] * alphaA[2];
	v[175] = 0.5e0*alphaA[0] * alphaA[2];
	v[166] = Power(alphaA[2], 2);
	v[774] = v[165] + v[166];
	v[188] = Power(alphaB[0], 2);
	v[186] = 0.5e0*alphaB[0] * alphaB[1];
	v[181] = Power(alphaB[1], 2);
	v[193] = 0.5e0*alphaB[1] * alphaB[2];
	v[191] = 0.5e0*alphaB[0] * alphaB[2];
	v[182] = Power(alphaB[2], 2);
	v[775] = v[181] + v[182];
	v[164] = 4e0 / (4e0 + v[172] + v[774]);
	v[167] = 1e0 - 0.5e0*v[164] * v[774];
	v[168] = v[164] * (-alphaA[2] + v[170]);
	v[169] = v[164] * (alphaA[1] + v[175]);
	v[171] = v[164] * (alphaA[2] + v[170]);
	v[173] = 1e0 - 0.5e0*v[164] * (v[166] + v[172]);
	v[174] = v[164] * (-alphaA[0] + v[177]);
	v[176] = v[164] * (-alphaA[1] + v[175]);
	v[178] = v[164] * (alphaA[0] + v[177]);
	v[179] = 1e0 - 0.5e0*v[164] * (v[165] + v[172]);
	v[180] = 4e0 / (4e0 + v[188] + v[775]);
	v[183] = 1e0 - 0.5e0*v[180] * v[775];
	v[184] = v[180] * (-alphaB[2] + v[186]);
	v[185] = v[180] * (alphaB[1] + v[191]);
	v[187] = v[180] * (alphaB[2] + v[186]);
	v[189] = 1e0 - 0.5e0*v[180] * (v[182] + v[188]);
	v[190] = v[180] * (-alphaB[0] + v[193]);
	v[192] = v[180] * (-alphaB[1] + v[191]);
	v[194] = v[180] * (alphaB[0] + v[193]);
	v[195] = 1e0 - 0.5e0*v[180] * (v[181] + v[188]);
	(Ob) = 0.5e0*(Power(uA[0] - uB[0] + GAp[0] * (QAi[0][0] * v[167] + QAi[1][0] * v[168] + QAi[2][0] * v[169]) + GAp[1] *
		(QAi[0][1] * v[167] + QAi[1][1] * v[168] + QAi[2][1] * v[169]) + GAp[2] * (QAi[0][2] * v[167] + QAi[1][2] * v[168]
			+ QAi[2][2] * v[169]) - GBp[0] * (QBi[0][0] * v[183] + QBi[1][0] * v[184] + QBi[2][0] * v[185]) - GBp[1] *
			(QBi[0][1] * v[183] + QBi[1][1] * v[184] + QBi[2][1] * v[185]) - GBp[2] * (QBi[0][2] * v[183] + QBi[1][2] * v[184]
				+ QBi[2][2] * v[185]) + xAi[0] - xBi[0], 2) + Power(uA[1] - uB[1] + GAp[0] * (QAi[0][0] * v[171] + QAi[1][0] * v[173]
					+ QAi[2][0] * v[174]) + GAp[1] * (QAi[0][1] * v[171] + QAi[1][1] * v[173] + QAi[2][1] * v[174]) + GAp[2] *
					(QAi[0][2] * v[171] + QAi[1][2] * v[173] + QAi[2][2] * v[174]) - GBp[0] * (QBi[0][0] * v[187] + QBi[1][0] * v[189]
						+ QBi[2][0] * v[190]) - GBp[1] * (QBi[0][1] * v[187] + QBi[1][1] * v[189] + QBi[2][1] * v[190]) - GBp[2] *
						(QBi[0][2] * v[187] + QBi[1][2] * v[189] + QBi[2][2] * v[190]) + xAi[1] - xBi[1], 2) + Power(uA[2] - uB[2] + GAp[0] *
					(QAi[0][0] * v[176] + QAi[1][0] * v[178] + QAi[2][0] * v[179]) + GAp[1] * (QAi[0][1] * v[176] + QAi[1][1] * v[178]
						+ QAi[2][1] * v[179]) + GAp[2] * (QAi[0][2] * v[176] + QAi[1][2] * v[178] + QAi[2][2] * v[179]) - GBp[0] *
						(QBi[0][0] * v[192] + QBi[1][0] * v[194] + QBi[2][0] * v[195]) - GBp[1] * (QBi[0][1] * v[192] + QBi[1][1] * v[194]
							+ QBi[2][1] * v[195]) - GBp[2] * (QBi[0][2] * v[192] + QBi[1][2] * v[194] + QBi[2][2] * v[195]) + xAi[2] - xBi[2], 2));

	return Ob;
}

//Calcula o Gradiente da função objetivo - Phase 1
void RigidNURBS_1_RigidNURBS_1::GradientPhase1(Matrix& mc, Matrix& mGra)
{
	//AceGen variables or pointers
	double v[2000];
	EvaluateNURBSDerivatives_p(mc);
	EvaluateNURBSDOFsVariables();

	double Gra[4];

	v[172] = Power(alphaA[0], 2);
	v[170] = 0.5e0*alphaA[0] * alphaA[1];
	v[165] = Power(alphaA[1], 2);
	v[177] = 0.5e0*alphaA[1] * alphaA[2];
	v[175] = 0.5e0*alphaA[0] * alphaA[2];
	v[166] = Power(alphaA[2], 2);
	v[774] = v[165] + v[166];
	v[188] = Power(alphaB[0], 2);
	v[186] = 0.5e0*alphaB[0] * alphaB[1];
	v[181] = Power(alphaB[1], 2);
	v[193] = 0.5e0*alphaB[1] * alphaB[2];
	v[191] = 0.5e0*alphaB[0] * alphaB[2];
	v[182] = Power(alphaB[2], 2);
	v[775] = v[181] + v[182];
	v[164] = 4e0 / (4e0 + v[172] + v[774]);
	v[167] = 1e0 - 0.5e0*v[164] * v[774];
	v[168] = v[164] * (-alphaA[2] + v[170]);
	v[169] = v[164] * (alphaA[1] + v[175]);
	v[171] = v[164] * (alphaA[2] + v[170]);
	v[173] = 1e0 - 0.5e0*v[164] * (v[166] + v[172]);
	v[174] = v[164] * (-alphaA[0] + v[177]);
	v[176] = v[164] * (-alphaA[1] + v[175]);
	v[178] = v[164] * (alphaA[0] + v[177]);
	v[179] = 1e0 - 0.5e0*v[164] * (v[165] + v[172]);
	v[180] = 4e0 / (4e0 + v[188] + v[775]);
	v[183] = 1e0 - 0.5e0*v[180] * v[775];
	v[184] = v[180] * (-alphaB[2] + v[186]);
	v[185] = v[180] * (alphaB[1] + v[191]);
	v[187] = v[180] * (alphaB[2] + v[186]);
	v[189] = 1e0 - 0.5e0*v[180] * (v[182] + v[188]);
	v[190] = v[180] * (-alphaB[0] + v[193]);
	v[192] = v[180] * (-alphaB[1] + v[191]);
	v[194] = v[180] * (alphaB[0] + v[193]);
	v[195] = 1e0 - 0.5e0*v[180] * (v[181] + v[188]);
	v[196] = QAi[0][0] * v[167] + QAi[1][0] * v[168] + QAi[2][0] * v[169];
	v[197] = QAi[0][1] * v[167] + QAi[1][1] * v[168] + QAi[2][1] * v[169];
	v[198] = QAi[0][2] * v[167] + QAi[1][2] * v[168] + QAi[2][2] * v[169];
	v[199] = QAi[0][0] * v[171] + QAi[1][0] * v[173] + QAi[2][0] * v[174];
	v[200] = QAi[0][1] * v[171] + QAi[1][1] * v[173] + QAi[2][1] * v[174];
	v[201] = QAi[0][2] * v[171] + QAi[1][2] * v[173] + QAi[2][2] * v[174];
	v[202] = QAi[0][0] * v[176] + QAi[1][0] * v[178] + QAi[2][0] * v[179];
	v[203] = QAi[0][1] * v[176] + QAi[1][1] * v[178] + QAi[2][1] * v[179];
	v[204] = QAi[0][2] * v[176] + QAi[1][2] * v[178] + QAi[2][2] * v[179];
	v[205] = QBi[0][0] * v[183] + QBi[1][0] * v[184] + QBi[2][0] * v[185];
	v[206] = QBi[0][1] * v[183] + QBi[1][1] * v[184] + QBi[2][1] * v[185];
	v[207] = QBi[0][2] * v[183] + QBi[1][2] * v[184] + QBi[2][2] * v[185];
	v[208] = QBi[0][0] * v[187] + QBi[1][0] * v[189] + QBi[2][0] * v[190];
	v[209] = QBi[0][1] * v[187] + QBi[1][1] * v[189] + QBi[2][1] * v[190];
	v[210] = QBi[0][2] * v[187] + QBi[1][2] * v[189] + QBi[2][2] * v[190];
	v[211] = QBi[0][0] * v[192] + QBi[1][0] * v[194] + QBi[2][0] * v[195];
	v[212] = QBi[0][1] * v[192] + QBi[1][1] * v[194] + QBi[2][1] * v[195];
	v[213] = QBi[0][2] * v[192] + QBi[1][2] * v[194] + QBi[2][2] * v[195];
	v[776] = 2e0*(uA[0] - uB[0] + GAp[0] * v[196] + GAp[1] * v[197] + GAp[2] * v[198] - GBp[0] * v[205] - GBp[1] * v[206]
		- GBp[2] * v[207] + xAi[0] - xBi[0]);
	v[777] = 2e0*(uA[1] - uB[1] + GAp[0] * v[199] + GAp[1] * v[200] + GAp[2] * v[201] - GBp[0] * v[208] - GBp[1] * v[209]
		- GBp[2] * v[210] + xAi[1] - xBi[1]);
	v[778] = 2e0*(uA[2] - uB[2] + GAp[0] * v[202] + GAp[1] * v[203] + GAp[2] * v[204] - GBp[0] * v[211] - GBp[1] * v[212]
		- GBp[2] * v[213] + xAi[2] - xBi[2]);
	Gra[0] = 0.5e0*((dGAp[0][0] * v[196] + dGAp[1][0] * v[197] + dGAp[2][0] * v[198])*v[776] + (dGAp[0][0] * v[199]
		+ dGAp[1][0] * v[200] + dGAp[2][0] * v[201])*v[777] + (dGAp[0][0] * v[202] + dGAp[1][0] * v[203] + dGAp[2][0] * v[204]
			)*v[778]);
	Gra[1] = 0.5e0*((dGAp[0][1] * v[196] + dGAp[1][1] * v[197] + dGAp[2][1] * v[198])*v[776] + (dGAp[0][1] * v[199]
		+ dGAp[1][1] * v[200] + dGAp[2][1] * v[201])*v[777] + (dGAp[0][1] * v[202] + dGAp[1][1] * v[203] + dGAp[2][1] * v[204]
			)*v[778]);
	Gra[2] = 0.5e0*((-(dGBp[0][0] * v[205]) - dGBp[1][0] * v[206] - dGBp[2][0] * v[207])*v[776] + (-
		(dGBp[0][0] * v[208]) - dGBp[1][0] * v[209] - dGBp[2][0] * v[210])*v[777] + (-(dGBp[0][0] * v[211])
			- dGBp[1][0] * v[212] - dGBp[2][0] * v[213])*v[778]);
	Gra[3] = 0.5e0*((-(dGBp[0][1] * v[205]) - dGBp[1][1] * v[206] - dGBp[2][1] * v[207])*v[776] + (-
		(dGBp[0][1] * v[208]) - dGBp[1][1] * v[209] - dGBp[2][1] * v[210])*v[777] + (-(dGBp[0][1] * v[211])
			- dGBp[1][1] * v[212] - dGBp[2][1] * v[213])*v[778]);

	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];
}

//Calcula a Hessiana da função objetivo - Phase 1
void RigidNURBS_1_RigidNURBS_1::HessianPhase1(Matrix& mc, Matrix& mHes)
{
	//AceGen variables or pointers
	double v[2000];
	EvaluateNURBSDerivatives_p(mc);
	EvaluateNURBSDOFsVariables();
	
	double Hes[4][4];

	v[172] = Power(alphaA[0], 2);
	v[170] = 0.5e0*alphaA[0] * alphaA[1];
	v[165] = Power(alphaA[1], 2);
	v[177] = 0.5e0*alphaA[1] * alphaA[2];
	v[175] = 0.5e0*alphaA[0] * alphaA[2];
	v[166] = Power(alphaA[2], 2);
	v[774] = v[165] + v[166];
	v[188] = Power(alphaB[0], 2);
	v[186] = 0.5e0*alphaB[0] * alphaB[1];
	v[181] = Power(alphaB[1], 2);
	v[193] = 0.5e0*alphaB[1] * alphaB[2];
	v[191] = 0.5e0*alphaB[0] * alphaB[2];
	v[182] = Power(alphaB[2], 2);
	v[775] = v[181] + v[182];
	v[164] = 4e0 / (4e0 + v[172] + v[774]);
	v[167] = 1e0 - 0.5e0*v[164] * v[774];
	v[168] = v[164] * (-alphaA[2] + v[170]);
	v[169] = v[164] * (alphaA[1] + v[175]);
	v[171] = v[164] * (alphaA[2] + v[170]);
	v[173] = 1e0 - 0.5e0*v[164] * (v[166] + v[172]);
	v[174] = v[164] * (-alphaA[0] + v[177]);
	v[176] = v[164] * (-alphaA[1] + v[175]);
	v[178] = v[164] * (alphaA[0] + v[177]);
	v[179] = 1e0 - 0.5e0*v[164] * (v[165] + v[172]);
	v[180] = 4e0 / (4e0 + v[188] + v[775]);
	v[183] = 1e0 - 0.5e0*v[180] * v[775];
	v[184] = v[180] * (-alphaB[2] + v[186]);
	v[185] = v[180] * (alphaB[1] + v[191]);
	v[187] = v[180] * (alphaB[2] + v[186]);
	v[189] = 1e0 - 0.5e0*v[180] * (v[182] + v[188]);
	v[190] = v[180] * (-alphaB[0] + v[193]);
	v[192] = v[180] * (-alphaB[1] + v[191]);
	v[194] = v[180] * (alphaB[0] + v[193]);
	v[195] = 1e0 - 0.5e0*v[180] * (v[181] + v[188]);
	v[196] = QAi[0][0] * v[167] + QAi[1][0] * v[168] + QAi[2][0] * v[169];
	v[197] = QAi[0][1] * v[167] + QAi[1][1] * v[168] + QAi[2][1] * v[169];
	v[198] = QAi[0][2] * v[167] + QAi[1][2] * v[168] + QAi[2][2] * v[169];
	v[228] = dGAp[0][1] * v[196] + dGAp[1][1] * v[197] + dGAp[2][1] * v[198];
	v[227] = dGAp[0][0] * v[196] + dGAp[1][0] * v[197] + dGAp[2][0] * v[198];
	v[199] = QAi[0][0] * v[171] + QAi[1][0] * v[173] + QAi[2][0] * v[174];
	v[200] = QAi[0][1] * v[171] + QAi[1][1] * v[173] + QAi[2][1] * v[174];
	v[201] = QAi[0][2] * v[171] + QAi[1][2] * v[173] + QAi[2][2] * v[174];
	v[230] = dGAp[0][1] * v[199] + dGAp[1][1] * v[200] + dGAp[2][1] * v[201];
	v[229] = dGAp[0][0] * v[199] + dGAp[1][0] * v[200] + dGAp[2][0] * v[201];
	v[202] = QAi[0][0] * v[176] + QAi[1][0] * v[178] + QAi[2][0] * v[179];
	v[203] = QAi[0][1] * v[176] + QAi[1][1] * v[178] + QAi[2][1] * v[179];
	v[204] = QAi[0][2] * v[176] + QAi[1][2] * v[178] + QAi[2][2] * v[179];
	v[232] = dGAp[0][1] * v[202] + dGAp[1][1] * v[203] + dGAp[2][1] * v[204];
	v[231] = dGAp[0][0] * v[202] + dGAp[1][0] * v[203] + dGAp[2][0] * v[204];
	v[282] = 2e0*(v[227] * v[228] + v[229] * v[230] + v[231] * v[232]);
	v[205] = QBi[0][0] * v[183] + QBi[1][0] * v[184] + QBi[2][0] * v[185];
	v[206] = QBi[0][1] * v[183] + QBi[1][1] * v[184] + QBi[2][1] * v[185];
	v[207] = QBi[0][2] * v[183] + QBi[1][2] * v[184] + QBi[2][2] * v[185];
	v[234] = dGBp[0][1] * v[205] + dGBp[1][1] * v[206] + dGBp[2][1] * v[207];
	v[233] = dGBp[0][0] * v[205] + dGBp[1][0] * v[206] + dGBp[2][0] * v[207];
	v[208] = QBi[0][0] * v[187] + QBi[1][0] * v[189] + QBi[2][0] * v[190];
	v[209] = QBi[0][1] * v[187] + QBi[1][1] * v[189] + QBi[2][1] * v[190];
	v[210] = QBi[0][2] * v[187] + QBi[1][2] * v[189] + QBi[2][2] * v[190];
	v[236] = dGBp[0][1] * v[208] + dGBp[1][1] * v[209] + dGBp[2][1] * v[210];
	v[235] = dGBp[0][0] * v[208] + dGBp[1][0] * v[209] + dGBp[2][0] * v[210];
	v[211] = QBi[0][0] * v[192] + QBi[1][0] * v[194] + QBi[2][0] * v[195];
	v[212] = QBi[0][1] * v[192] + QBi[1][1] * v[194] + QBi[2][1] * v[195];
	v[213] = QBi[0][2] * v[192] + QBi[1][2] * v[194] + QBi[2][2] * v[195];
	v[238] = dGBp[0][1] * v[211] + dGBp[1][1] * v[212] + dGBp[2][1] * v[213];
	v[271] = -1e0*(v[228] * v[234] + v[230] * v[236] + v[232] * v[238]);
	v[270] = -1e0*(v[227] * v[234] + v[229] * v[236] + v[231] * v[238]);
	v[237] = dGBp[0][0] * v[211] + dGBp[1][0] * v[212] + dGBp[2][0] * v[213];
	v[277] = 2e0*(v[233] * v[234] + v[235] * v[236] + v[237] * v[238]);
	v[275] = -1e0*(v[228] * v[233] + v[230] * v[235] + v[232] * v[237]);
	v[274] = -1e0*(v[227] * v[233] + v[229] * v[235] + v[231] * v[237]);
	v[776] = 2e0*(uA[0] - uB[0] + GAp[0] * v[196] + GAp[1] * v[197] + GAp[2] * v[198] - GBp[0] * v[205] - GBp[1] * v[206]
		- GBp[2] * v[207] + xAi[0] - xBi[0]);
	v[777] = 2e0*(uA[1] - uB[1] + GAp[0] * v[199] + GAp[1] * v[200] + GAp[2] * v[201] - GBp[0] * v[208] - GBp[1] * v[209]
		- GBp[2] * v[210] + xAi[1] - xBi[1]);
	v[778] = 2e0*(uA[2] - uB[2] + GAp[0] * v[202] + GAp[1] * v[203] + GAp[2] * v[204] - GBp[0] * v[211] - GBp[1] * v[212]
		- GBp[2] * v[213] + xAi[2] - xBi[2]);
	Hes[0][0] = 0.5e0*(2e0*(v[227] * v[227]) + 2e0*(v[229] * v[229]) + 2e0*(v[231] * v[231]) + (ddGAp[0][0][0] * v[196]
		+ ddGAp[1][0][0] * v[197] + ddGAp[2][0][0] * v[198])*v[776] + (ddGAp[0][0][0] * v[199] + ddGAp[1][0][0] * v[200]
			+ ddGAp[2][0][0] * v[201])*v[777] + (ddGAp[0][0][0] * v[202] + ddGAp[1][0][0] * v[203] + ddGAp[2][0][0] * v[204]
				)*v[778]);
	Hes[0][1] = 0.5e0*(v[282] + (ddGAp[0][0][1] * v[196] + ddGAp[1][0][1] * v[197] + ddGAp[2][0][1] * v[198])*v[776] +
		(ddGAp[0][0][1] * v[199] + ddGAp[1][0][1] * v[200] + ddGAp[2][0][1] * v[201])*v[777] + (ddGAp[0][0][1] * v[202]
			+ ddGAp[1][0][1] * v[203] + ddGAp[2][0][1] * v[204])*v[778]);
	Hes[0][2] = v[274];
	Hes[0][3] = v[270];
	Hes[1][0] = 0.5e0*(v[282] + (ddGAp[0][1][0] * v[196] + ddGAp[1][1][0] * v[197] + ddGAp[2][1][0] * v[198])*v[776] +
		(ddGAp[0][1][0] * v[199] + ddGAp[1][1][0] * v[200] + ddGAp[2][1][0] * v[201])*v[777] + (ddGAp[0][1][0] * v[202]
			+ ddGAp[1][1][0] * v[203] + ddGAp[2][1][0] * v[204])*v[778]);
	Hes[1][1] = 0.5e0*(2e0*(v[228] * v[228]) + 2e0*(v[230] * v[230]) + 2e0*(v[232] * v[232]) + (ddGAp[0][1][1] * v[196]
		+ ddGAp[1][1][1] * v[197] + ddGAp[2][1][1] * v[198])*v[776] + (ddGAp[0][1][1] * v[199] + ddGAp[1][1][1] * v[200]
			+ ddGAp[2][1][1] * v[201])*v[777] + (ddGAp[0][1][1] * v[202] + ddGAp[1][1][1] * v[203] + ddGAp[2][1][1] * v[204]
				)*v[778]);
	Hes[1][2] = v[275];
	Hes[1][3] = v[271];
	Hes[2][0] = v[274];
	Hes[2][1] = v[275];
	Hes[2][2] = 0.5e0*(2e0*(v[233] * v[233]) + 2e0*(v[235] * v[235]) + 2e0*(v[237] * v[237]) + (-
		(ddGBp[0][0][0] * v[205]) - ddGBp[1][0][0] * v[206] - ddGBp[2][0][0] * v[207])*v[776] + (-
		(ddGBp[0][0][0] * v[208]) - ddGBp[1][0][0] * v[209] - ddGBp[2][0][0] * v[210])*v[777] + (-
			(ddGBp[0][0][0] * v[211]) - ddGBp[1][0][0] * v[212] - ddGBp[2][0][0] * v[213])*v[778]);
	Hes[2][3] = 0.5e0*(v[277] + (-(ddGBp[0][0][1] * v[205]) - ddGBp[1][0][1] * v[206] - ddGBp[2][0][1] * v[207]
		)*v[776] + (-(ddGBp[0][0][1] * v[208]) - ddGBp[1][0][1] * v[209] - ddGBp[2][0][1] * v[210])*v[777] + (-
		(ddGBp[0][0][1] * v[211]) - ddGBp[1][0][1] * v[212] - ddGBp[2][0][1] * v[213])*v[778]);
	Hes[3][0] = v[270];
	Hes[3][1] = v[271];
	Hes[3][2] = 0.5e0*(v[277] + (-(ddGBp[0][1][0] * v[205]) - ddGBp[1][1][0] * v[206] - ddGBp[2][1][0] * v[207]
		)*v[776] + (-(ddGBp[0][1][0] * v[208]) - ddGBp[1][1][0] * v[209] - ddGBp[2][1][0] * v[210])*v[777] + (-
		(ddGBp[0][1][0] * v[211]) - ddGBp[1][1][0] * v[212] - ddGBp[2][1][0] * v[213])*v[778]);
	Hes[3][3] = 0.5e0*(2e0*(v[234] * v[234]) + 2e0*(v[236] * v[236]) + 2e0*(v[238] * v[238]) + (-
		(ddGBp[0][1][1] * v[205]) - ddGBp[1][1][1] * v[206] - ddGBp[2][1][1] * v[207])*v[776] + (-
		(ddGBp[0][1][1] * v[208]) - ddGBp[1][1][1] * v[209] - ddGBp[2][1][1] * v[210])*v[777] + (-
			(ddGBp[0][1][1] * v[211]) - ddGBp[1][1][1] * v[212] - ddGBp[2][1][1] * v[213])*v[778]);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mHes(i, j) = Hes[i][j];
}

//Calcula e rotorna o gap (com sinal)
double RigidNURBS_1_RigidNURBS_1::Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	//AceGen variables or pointers
	double v[2000];
	EvaluateNURBSDerivatives_p(mc);
	EvaluateNURBSDOFsVariables();
	double Gap;
	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	int b243;
	v[161] = Power(alphaA[0], 2);
	v[159] = 0.5e0*alphaA[0] * alphaA[1];
	v[154] = Power(alphaA[1], 2);
	v[166] = 0.5e0*alphaA[1] * alphaA[2];
	v[164] = 0.5e0*alphaA[0] * alphaA[2];
	v[155] = Power(alphaA[2], 2);
	v[252] = v[154] + v[155];
	v[192] = Power(alphaB[0], 2);
	v[190] = 0.5e0*alphaB[0] * alphaB[1];
	v[185] = Power(alphaB[1], 2);
	v[197] = 0.5e0*alphaB[1] * alphaB[2];
	v[195] = 0.5e0*alphaB[0] * alphaB[2];
	v[186] = Power(alphaB[2], 2);
	v[253] = v[185] + v[186];
	v[153] = 4e0 / (4e0 + v[161] + v[252]);
	v[156] = 1e0 - 0.5e0*v[153] * v[252];
	v[157] = v[153] * (-alphaA[2] + v[159]);
	v[158] = v[153] * (alphaA[1] + v[164]);
	v[160] = v[153] * (alphaA[2] + v[159]);
	v[162] = 1e0 - 0.5e0*v[153] * (v[155] + v[161]);
	v[163] = v[153] * (-alphaA[0] + v[166]);
	v[165] = v[153] * (-alphaA[1] + v[164]);
	v[167] = v[153] * (alphaA[0] + v[166]);
	v[168] = 1e0 - 0.5e0*v[153] * (v[154] + v[161]);
	v[169] = QAi[0][0] * v[156] + QAi[1][0] * v[157] + QAi[2][0] * v[158];
	v[170] = QAi[0][1] * v[156] + QAi[1][1] * v[157] + QAi[2][1] * v[158];
	v[171] = QAi[0][2] * v[156] + QAi[1][2] * v[157] + QAi[2][2] * v[158];
	v[218] = dGAp[0][1] * v[169] + dGAp[1][1] * v[170] + dGAp[2][1] * v[171];
	v[215] = dGAp[0][0] * v[169] + dGAp[1][0] * v[170] + dGAp[2][0] * v[171];
	v[172] = QAi[0][0] * v[160] + QAi[1][0] * v[162] + QAi[2][0] * v[163];
	v[173] = QAi[0][1] * v[160] + QAi[1][1] * v[162] + QAi[2][1] * v[163];
	v[174] = QAi[0][2] * v[160] + QAi[1][2] * v[162] + QAi[2][2] * v[163];
	v[219] = dGAp[0][1] * v[172] + dGAp[1][1] * v[173] + dGAp[2][1] * v[174];
	v[216] = dGAp[0][0] * v[172] + dGAp[1][0] * v[173] + dGAp[2][0] * v[174];
	v[175] = QAi[0][0] * v[165] + QAi[1][0] * v[167] + QAi[2][0] * v[168];
	v[176] = QAi[0][1] * v[165] + QAi[1][1] * v[167] + QAi[2][1] * v[168];
	v[177] = QAi[0][2] * v[165] + QAi[1][2] * v[167] + QAi[2][2] * v[168];
	v[220] = dGAp[0][1] * v[175] + dGAp[1][1] * v[176] + dGAp[2][1] * v[177];
	v[217] = dGAp[0][0] * v[175] + dGAp[1][0] * v[176] + dGAp[2][0] * v[177];
	v[184] = 4e0 / (4e0 + v[192] + v[253]);
	v[187] = 1e0 - 0.5e0*v[184] * v[253];
	v[188] = v[184] * (-alphaB[2] + v[190]);
	v[189] = v[184] * (alphaB[1] + v[195]);
	v[191] = v[184] * (alphaB[2] + v[190]);
	v[193] = 1e0 - 0.5e0*v[184] * (v[186] + v[192]);
	v[194] = v[184] * (-alphaB[0] + v[197]);
	v[196] = v[184] * (-alphaB[1] + v[195]);
	v[198] = v[184] * (alphaB[0] + v[197]);
	v[199] = 1e0 - 0.5e0*v[184] * (v[185] + v[192]);
	v[200] = QBi[0][0] * v[187] + QBi[1][0] * v[188] + QBi[2][0] * v[189];
	v[201] = QBi[0][1] * v[187] + QBi[1][1] * v[188] + QBi[2][1] * v[189];
	v[202] = QBi[0][2] * v[187] + QBi[1][2] * v[188] + QBi[2][2] * v[189];
	v[224] = dGBp[0][1] * v[200] + dGBp[1][1] * v[201] + dGBp[2][1] * v[202];
	v[221] = dGBp[0][0] * v[200] + dGBp[1][0] * v[201] + dGBp[2][0] * v[202];
	v[203] = QBi[0][0] * v[191] + QBi[1][0] * v[193] + QBi[2][0] * v[194];
	v[204] = QBi[0][1] * v[191] + QBi[1][1] * v[193] + QBi[2][1] * v[194];
	v[205] = QBi[0][2] * v[191] + QBi[1][2] * v[193] + QBi[2][2] * v[194];
	v[225] = dGBp[0][1] * v[203] + dGBp[1][1] * v[204] + dGBp[2][1] * v[205];
	v[222] = dGBp[0][0] * v[203] + dGBp[1][0] * v[204] + dGBp[2][0] * v[205];
	v[206] = QBi[0][0] * v[196] + QBi[1][0] * v[198] + QBi[2][0] * v[199];
	v[207] = QBi[0][1] * v[196] + QBi[1][1] * v[198] + QBi[2][1] * v[199];
	v[208] = QBi[0][2] * v[196] + QBi[1][2] * v[198] + QBi[2][2] * v[199];
	v[226] = dGBp[0][1] * v[206] + dGBp[1][1] * v[207] + dGBp[2][1] * v[208];
	v[223] = dGBp[0][0] * v[206] + dGBp[1][0] * v[207] + dGBp[2][0] * v[208];
	v[248] = uA[0] - uB[0] + GAp[0] * v[169] + GAp[1] * v[170] + GAp[2] * v[171] - GBp[0] * v[200] - GBp[1] * v[201]
		- GBp[2] * v[202] + xAi[0] - xBi[0];
	v[249] = uA[1] - uB[1] + GAp[0] * v[172] + GAp[1] * v[173] + GAp[2] * v[174] - GBp[0] * v[203] - GBp[1] * v[204]
		- GBp[2] * v[205] + xAi[1] - xBi[1];
	v[250] = uA[2] - uB[2] + GAp[0] * v[175] + GAp[1] * v[176] + GAp[2] * v[177] - GBp[0] * v[206] - GBp[1] * v[207]
		- GBp[2] * v[208] + xAi[2] - xBi[2];
	v[229] = -(v[217] * v[219]) + v[216] * v[220];
	v[236] = -(v[223] * v[225]) + v[222] * v[226];
	if ((*fixnormal)) {
		v[247] = ((-normalA[0] + normalB[0])*v[248] + (-normalA[1] + normalB[1])*v[249] + (-normalA[2] + normalB[2]
			)*v[250]) / 2e0;
	}
	else {
		v[241] = -(v[222] * v[224]) + v[221] * v[225];
		v[239] = v[223] * v[224] - v[221] * v[226];
		v[234] = -(v[216] * v[218]) + v[215] * v[219];
		v[232] = v[217] * v[218] - v[215] * v[220];
		v[247] = (-((((*invertnormalA) ? -1 : 1)*(v[229] * v[248] + v[232] * v[249] + v[234] * v[250])) / sqrt(
			(v[229] * v[229]) + (v[232] * v[232]) + (v[234] * v[234]))) + (((*invertnormalB) ? -1 : 1)*(v[236] * v[248]
				+ v[239] * v[249] + v[241] * v[250])) / sqrt((v[236] * v[236]) + (v[239] * v[239]) + (v[241] * v[241]))) / 2e0;
	};
	(Gap) = v[247];
	

	return Gap;
}

//Calcula o Gradiente do gap
void RigidNURBS_1_RigidNURBS_1::GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	//AceGen variables or pointers
	double v[2000];
	EvaluateNURBSDerivatives_p(mc);
	EvaluateNURBSDOFsVariables();
	double Gra[4];
	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();
	
	int i252, i312, i313, b243, b254;
	v[161] = Power(alphaA[0], 2);
	v[159] = 0.5e0*alphaA[0] * alphaA[1];
	v[154] = Power(alphaA[1], 2);
	v[166] = 0.5e0*alphaA[1] * alphaA[2];
	v[164] = 0.5e0*alphaA[0] * alphaA[2];
	v[155] = Power(alphaA[2], 2);
	v[305] = v[154] + v[155];
	v[192] = Power(alphaB[0], 2);
	v[190] = 0.5e0*alphaB[0] * alphaB[1];
	v[185] = Power(alphaB[1], 2);
	v[197] = 0.5e0*alphaB[1] * alphaB[2];
	v[195] = 0.5e0*alphaB[0] * alphaB[2];
	v[186] = Power(alphaB[2], 2);
	v[306] = v[185] + v[186];
	v[153] = 4e0 / (4e0 + v[161] + v[305]);
	v[156] = 1e0 - 0.5e0*v[153] * v[305];
	v[157] = v[153] * (-alphaA[2] + v[159]);
	v[158] = v[153] * (alphaA[1] + v[164]);
	v[160] = v[153] * (alphaA[2] + v[159]);
	v[162] = 1e0 - 0.5e0*v[153] * (v[155] + v[161]);
	v[163] = v[153] * (-alphaA[0] + v[166]);
	v[165] = v[153] * (-alphaA[1] + v[164]);
	v[167] = v[153] * (alphaA[0] + v[166]);
	v[168] = 1e0 - 0.5e0*v[153] * (v[154] + v[161]);
	v[169] = QAi[0][0] * v[156] + QAi[1][0] * v[157] + QAi[2][0] * v[158];
	v[170] = QAi[0][1] * v[156] + QAi[1][1] * v[157] + QAi[2][1] * v[158];
	v[171] = QAi[0][2] * v[156] + QAi[1][2] * v[157] + QAi[2][2] * v[158];
	v[218] = dGAp[0][1] * v[169] + dGAp[1][1] * v[170] + dGAp[2][1] * v[171];
	v[215] = dGAp[0][0] * v[169] + dGAp[1][0] * v[170] + dGAp[2][0] * v[171];
	v[172] = QAi[0][0] * v[160] + QAi[1][0] * v[162] + QAi[2][0] * v[163];
	v[173] = QAi[0][1] * v[160] + QAi[1][1] * v[162] + QAi[2][1] * v[163];
	v[174] = QAi[0][2] * v[160] + QAi[1][2] * v[162] + QAi[2][2] * v[163];
	v[219] = dGAp[0][1] * v[172] + dGAp[1][1] * v[173] + dGAp[2][1] * v[174];
	v[216] = dGAp[0][0] * v[172] + dGAp[1][0] * v[173] + dGAp[2][0] * v[174];
	v[234] = -(v[216] * v[218]) + v[215] * v[219];
	v[175] = QAi[0][0] * v[165] + QAi[1][0] * v[167] + QAi[2][0] * v[168];
	v[176] = QAi[0][1] * v[165] + QAi[1][1] * v[167] + QAi[2][1] * v[168];
	v[177] = QAi[0][2] * v[165] + QAi[1][2] * v[167] + QAi[2][2] * v[168];
	v[220] = dGAp[0][1] * v[175] + dGAp[1][1] * v[176] + dGAp[2][1] * v[177];
	v[217] = dGAp[0][0] * v[175] + dGAp[1][0] * v[176] + dGAp[2][0] * v[177];
	v[232] = v[217] * v[218] - v[215] * v[220];
	v[184] = 4e0 / (4e0 + v[192] + v[306]);
	v[187] = 1e0 - 0.5e0*v[184] * v[306];
	v[188] = v[184] * (-alphaB[2] + v[190]);
	v[189] = v[184] * (alphaB[1] + v[195]);
	v[191] = v[184] * (alphaB[2] + v[190]);
	v[193] = 1e0 - 0.5e0*v[184] * (v[186] + v[192]);
	v[194] = v[184] * (-alphaB[0] + v[197]);
	v[196] = v[184] * (-alphaB[1] + v[195]);
	v[198] = v[184] * (alphaB[0] + v[197]);
	v[199] = 1e0 - 0.5e0*v[184] * (v[185] + v[192]);
	v[200] = QBi[0][0] * v[187] + QBi[1][0] * v[188] + QBi[2][0] * v[189];
	v[201] = QBi[0][1] * v[187] + QBi[1][1] * v[188] + QBi[2][1] * v[189];
	v[202] = QBi[0][2] * v[187] + QBi[1][2] * v[188] + QBi[2][2] * v[189];
	v[224] = dGBp[0][1] * v[200] + dGBp[1][1] * v[201] + dGBp[2][1] * v[202];
	v[221] = dGBp[0][0] * v[200] + dGBp[1][0] * v[201] + dGBp[2][0] * v[202];
	v[203] = QBi[0][0] * v[191] + QBi[1][0] * v[193] + QBi[2][0] * v[194];
	v[204] = QBi[0][1] * v[191] + QBi[1][1] * v[193] + QBi[2][1] * v[194];
	v[205] = QBi[0][2] * v[191] + QBi[1][2] * v[193] + QBi[2][2] * v[194];
	v[225] = dGBp[0][1] * v[203] + dGBp[1][1] * v[204] + dGBp[2][1] * v[205];
	v[222] = dGBp[0][0] * v[203] + dGBp[1][0] * v[204] + dGBp[2][0] * v[205];
	v[241] = -(v[222] * v[224]) + v[221] * v[225];
	v[206] = QBi[0][0] * v[196] + QBi[1][0] * v[198] + QBi[2][0] * v[199];
	v[207] = QBi[0][1] * v[196] + QBi[1][1] * v[198] + QBi[2][1] * v[199];
	v[208] = QBi[0][2] * v[196] + QBi[1][2] * v[198] + QBi[2][2] * v[199];
	v[226] = dGBp[0][1] * v[206] + dGBp[1][1] * v[207] + dGBp[2][1] * v[208];
	v[223] = dGBp[0][0] * v[206] + dGBp[1][0] * v[207] + dGBp[2][0] * v[208];
	v[239] = v[223] * v[224] - v[221] * v[226];
	v[248] = uA[0] - uB[0] + GAp[0] * v[169] + GAp[1] * v[170] + GAp[2] * v[171] - GBp[0] * v[200] - GBp[1] * v[201]
		- GBp[2] * v[202] + xAi[0] - xBi[0];
	v[249] = uA[1] - uB[1] + GAp[0] * v[172] + GAp[1] * v[173] + GAp[2] * v[174] - GBp[0] * v[203] - GBp[1] * v[204]
		- GBp[2] * v[205] + xAi[1] - xBi[1];
	v[250] = uA[2] - uB[2] + GAp[0] * v[175] + GAp[1] * v[176] + GAp[2] * v[177] - GBp[0] * v[206] - GBp[1] * v[207]
		- GBp[2] * v[208] + xAi[2] - xBi[2];
	v[229] = -(v[217] * v[219]) + v[216] * v[220];
	v[270] = (v[229] * v[229]) + (v[232] * v[232]) + (v[234] * v[234]);
	v[231] = 1e0 / sqrt(v[270]);
	v[236] = -(v[223] * v[225]) + v[222] * v[226];
	v[268] = (v[236] * v[236]) + (v[239] * v[239]) + (v[241] * v[241]);
	v[238] = 1e0 / sqrt(v[268]);
	if ((*fixnormal)) {
		v[311] = -normalA[0] + normalB[0];
		v[310] = -normalA[1] + normalB[1];
		v[309] = -normalA[2] + normalB[2];
		v[247] = (v[250] * v[309] + v[249] * v[310] + v[248] * v[311]) / 2e0;
	}
	else {
		i313 = ((*invertnormalA) ? -1 : 1);
		i312 = ((*invertnormalB) ? -1 : 1);
		v[308] = i313 * v[231];
		v[307] = i312 * v[238];
		v[314] = v[241] * v[307] - v[234] * v[308];
		v[315] = v[239] * v[307] - v[232] * v[308];
		v[316] = v[236] * v[307] - v[229] * v[308];
		v[247] = (v[250] * v[314] + v[249] * v[315] + v[248] * v[316]) / 2e0;
	};
	b254 = (*fixnormal);
	if (b254) {
		v[255] = 0e0;
		v[256] = 0e0;
		v[257] = 0e0;
		v[258] = 0e0;
		v[259] = 0e0;
		v[260] = 0e0;
		v[261] = 0e0;
		v[262] = 0e0;
		v[263] = v[309] / 2e0;
		v[264] = v[310] / 2e0;
		v[265] = v[311] / 2e0;
	}
	else {
		v[267] = -v[308] / 2e0;
		v[266] = v[307] / 2e0;
		v[260] = (i312*(v[236] * v[248] + v[239] * v[249] + v[241] * v[250])) / 2e0;
		v[259] = v[248] * v[266];
		v[256] = -(i313*(v[229] * v[248] + v[232] * v[249] + v[234] * v[250])) / 2e0;
		v[255] = v[248] * v[267];
		v[263] = v[314] / 2e0;
		v[264] = v[315] / 2e0;
		v[265] = v[316] / 2e0;
		v[261] = v[249] * v[266];
		v[262] = v[250] * v[266];
		v[257] = v[249] * v[267];
		v[258] = v[250] * v[267];
	};
	v[271] = -((v[231] * v[256]) / v[270]);
	v[269] = -((v[238] * v[260]) / v[268]);
	v[259] = v[259] + v[236] * v[269];
	v[261] = v[261] + v[239] * v[269];
	v[262] = v[262] + v[241] * v[269];
	v[255] = v[255] + v[229] * v[271];
	v[257] = v[257] + v[232] * v[271];
	v[258] = v[258] + v[234] * v[271];
	v[272] = -(v[208] * v[263]) - v[205] * v[264] - v[202] * v[265];
	v[273] = -(v[207] * v[263]) - v[204] * v[264] - v[201] * v[265];
	v[274] = -(v[206] * v[263]) - v[203] * v[264] - v[200] * v[265];
	v[275] = v[177] * v[263] + v[174] * v[264] + v[171] * v[265];
	v[276] = v[176] * v[263] + v[173] * v[264] + v[170] * v[265];
	v[277] = v[175] * v[263] + v[172] * v[264] + v[169] * v[265];
	v[278] = -(v[225] * v[259]) + v[224] * v[261];
	v[279] = v[222] * v[259] - v[221] * v[261];
	v[280] = v[226] * v[259] - v[224] * v[262];
	v[281] = -(v[223] * v[259]) + v[221] * v[262];
	v[282] = -(v[226] * v[261]) + v[225] * v[262];
	v[283] = v[223] * v[261] - v[222] * v[262];
	v[284] = v[208] * v[278] + v[205] * v[280] + v[202] * v[282];
	v[285] = v[207] * v[278] + v[204] * v[280] + v[201] * v[282];
	v[286] = v[206] * v[278] + v[203] * v[280] + v[200] * v[282];
	v[287] = v[208] * v[279] + v[205] * v[281] + v[202] * v[283];
	v[288] = v[207] * v[279] + v[204] * v[281] + v[201] * v[283];
	v[289] = v[206] * v[279] + v[203] * v[281] + v[200] * v[283];
	v[290] = -(v[219] * v[255]) + v[218] * v[257];
	v[291] = v[216] * v[255] - v[215] * v[257];
	v[292] = v[220] * v[255] - v[218] * v[258];
	v[293] = -(v[217] * v[255]) + v[215] * v[258];
	v[294] = -(v[220] * v[257]) + v[219] * v[258];
	v[295] = v[217] * v[257] - v[216] * v[258];
	v[296] = v[177] * v[290] + v[174] * v[292] + v[171] * v[294];
	v[297] = v[176] * v[290] + v[173] * v[292] + v[170] * v[294];
	v[298] = v[175] * v[290] + v[172] * v[292] + v[169] * v[294];
	v[299] = v[177] * v[291] + v[174] * v[293] + v[171] * v[295];
	v[300] = v[176] * v[291] + v[173] * v[293] + v[170] * v[295];
	v[301] = v[175] * v[291] + v[172] * v[293] + v[169] * v[295];
	v[365] = dGAp[2][0] * v[275] + dGAp[1][0] * v[276] + dGAp[0][0] * v[277] + ddGAp[2][0][0] * v[296]
		+ ddGAp[1][0][0] * v[297] + ddGAp[0][0][0] * v[298] + ddGAp[2][1][0] * v[299] + ddGAp[1][1][0] * v[300]
		+ ddGAp[0][1][0] * v[301];
	v[366] = dGAp[2][1] * v[275] + dGAp[1][1] * v[276] + dGAp[0][1] * v[277] + ddGAp[2][0][1] * v[296]
		+ ddGAp[1][0][1] * v[297] + ddGAp[0][0][1] * v[298] + ddGAp[2][1][1] * v[299] + ddGAp[1][1][1] * v[300]
		+ ddGAp[0][1][1] * v[301];
	v[367] = dGBp[2][0] * v[272] + dGBp[1][0] * v[273] + dGBp[0][0] * v[274] + ddGBp[2][0][0] * v[284]
		+ ddGBp[1][0][0] * v[285] + ddGBp[0][0][0] * v[286] + ddGBp[2][1][0] * v[287] + ddGBp[1][1][0] * v[288]
		+ ddGBp[0][1][0] * v[289];
	v[368] = dGBp[2][1] * v[272] + dGBp[1][1] * v[273] + dGBp[0][1] * v[274] + ddGBp[2][0][1] * v[284]
		+ ddGBp[1][0][1] * v[285] + ddGBp[0][0][1] * v[286] + ddGBp[2][1][1] * v[287] + ddGBp[1][1][1] * v[288]
		+ ddGBp[0][1][1] * v[289];
	
	for (i252 = 1; i252 <= 4; i252++) {
		Gra[i252 - 1] = v[364 + i252];
	};/* end for */

	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];
}

//Calcula a Hessiana do gap
void RigidNURBS_1_RigidNURBS_1::HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	//AceGen variables or pointers
	double v[2000];
	EvaluateNURBSDerivatives_p(mc);
	EvaluateNURBSDOFsVariables();
	double Hes[4][4];
	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();
	
	int i252, i305, i425, i426, i427, i428, b243, b254, b365, b411;
	v[608] = 0e0;
	v[609] = 0e0;
	v[610] = dGBp[2][0];
	v[611] = dGBp[2][1];
	v[604] = 0e0;
	v[605] = 0e0;
	v[606] = dGBp[1][0];
	v[607] = dGBp[1][1];
	v[600] = 0e0;
	v[601] = 0e0;
	v[602] = dGBp[0][0];
	v[603] = dGBp[0][1];
	v[596] = dGAp[2][0];
	v[597] = dGAp[2][1];
	v[598] = 0e0;
	v[599] = 0e0;
	v[592] = dGAp[1][0];
	v[593] = dGAp[1][1];
	v[594] = 0e0;
	v[595] = 0e0;
	v[588] = dGAp[0][0];
	v[589] = dGAp[0][1];
	v[590] = 0e0;
	v[591] = 0e0;
	v[584] = 0e0;
	v[585] = 0e0;
	v[586] = ddGBp[2][0][0];
	v[587] = ddGBp[2][0][1];
	v[580] = 0e0;
	v[581] = 0e0;
	v[582] = ddGBp[1][0][0];
	v[583] = ddGBp[1][0][1];
	v[576] = 0e0;
	v[577] = 0e0;
	v[578] = ddGBp[0][0][0];
	v[579] = ddGBp[0][0][1];
	v[572] = 0e0;
	v[573] = 0e0;
	v[574] = ddGBp[2][1][0];
	v[575] = ddGBp[2][1][1];
	v[568] = 0e0;
	v[569] = 0e0;
	v[570] = ddGBp[1][1][0];
	v[571] = ddGBp[1][1][1];
	v[564] = 0e0;
	v[565] = 0e0;
	v[566] = ddGBp[0][1][0];
	v[567] = ddGBp[0][1][1];
	v[560] = ddGAp[2][0][0];
	v[561] = ddGAp[2][0][1];
	v[562] = 0e0;
	v[563] = 0e0;
	v[556] = ddGAp[1][0][0];
	v[557] = ddGAp[1][0][1];
	v[558] = 0e0;
	v[559] = 0e0;
	v[552] = ddGAp[0][0][0];
	v[553] = ddGAp[0][0][1];
	v[554] = 0e0;
	v[555] = 0e0;
	v[548] = ddGAp[2][1][0];
	v[549] = ddGAp[2][1][1];
	v[550] = 0e0;
	v[551] = 0e0;
	v[544] = ddGAp[1][1][0];
	v[545] = ddGAp[1][1][1];
	v[546] = 0e0;
	v[547] = 0e0;
	v[540] = ddGAp[0][1][0];
	v[541] = ddGAp[0][1][1];
	v[542] = 0e0;
	v[543] = 0e0;
	v[161] = Power(alphaA[0], 2);
	v[159] = 0.5e0*alphaA[0] * alphaA[1];
	v[154] = Power(alphaA[1], 2);
	v[166] = 0.5e0*alphaA[1] * alphaA[2];
	v[164] = 0.5e0*alphaA[0] * alphaA[2];
	v[155] = Power(alphaA[2], 2);
	v[413] = v[154] + v[155];
	v[192] = Power(alphaB[0], 2);
	v[190] = 0.5e0*alphaB[0] * alphaB[1];
	v[185] = Power(alphaB[1], 2);
	v[197] = 0.5e0*alphaB[1] * alphaB[2];
	v[195] = 0.5e0*alphaB[0] * alphaB[2];
	v[186] = Power(alphaB[2], 2);
	v[414] = v[185] + v[186];
	v[153] = 4e0 / (4e0 + v[161] + v[413]);
	v[156] = 1e0 - 0.5e0*v[153] * v[413];
	v[157] = v[153] * (-alphaA[2] + v[159]);
	v[158] = v[153] * (alphaA[1] + v[164]);
	v[160] = v[153] * (alphaA[2] + v[159]);
	v[162] = 1e0 - 0.5e0*v[153] * (v[155] + v[161]);
	v[163] = v[153] * (-alphaA[0] + v[166]);
	v[165] = v[153] * (-alphaA[1] + v[164]);
	v[167] = v[153] * (alphaA[0] + v[166]);
	v[168] = 1e0 - 0.5e0*v[153] * (v[154] + v[161]);
	v[169] = QAi[0][0] * v[156] + QAi[1][0] * v[157] + QAi[2][0] * v[158];
	v[170] = QAi[0][1] * v[156] + QAi[1][1] * v[157] + QAi[2][1] * v[158];
	v[171] = QAi[0][2] * v[156] + QAi[1][2] * v[157] + QAi[2][2] * v[158];
	v[218] = dGAp[0][1] * v[169] + dGAp[1][1] * v[170] + dGAp[2][1] * v[171];
	v[215] = dGAp[0][0] * v[169] + dGAp[1][0] * v[170] + dGAp[2][0] * v[171];
	v[172] = QAi[0][0] * v[160] + QAi[1][0] * v[162] + QAi[2][0] * v[163];
	v[173] = QAi[0][1] * v[160] + QAi[1][1] * v[162] + QAi[2][1] * v[163];
	v[174] = QAi[0][2] * v[160] + QAi[1][2] * v[162] + QAi[2][2] * v[163];
	v[219] = dGAp[0][1] * v[172] + dGAp[1][1] * v[173] + dGAp[2][1] * v[174];
	v[216] = dGAp[0][0] * v[172] + dGAp[1][0] * v[173] + dGAp[2][0] * v[174];
	v[234] = -(v[216] * v[218]) + v[215] * v[219];
	v[175] = QAi[0][0] * v[165] + QAi[1][0] * v[167] + QAi[2][0] * v[168];
	v[176] = QAi[0][1] * v[165] + QAi[1][1] * v[167] + QAi[2][1] * v[168];
	v[177] = QAi[0][2] * v[165] + QAi[1][2] * v[167] + QAi[2][2] * v[168];
	v[220] = dGAp[0][1] * v[175] + dGAp[1][1] * v[176] + dGAp[2][1] * v[177];
	v[217] = dGAp[0][0] * v[175] + dGAp[1][0] * v[176] + dGAp[2][0] * v[177];
	v[232] = v[217] * v[218] - v[215] * v[220];
	v[184] = 4e0 / (4e0 + v[192] + v[414]);
	v[187] = 1e0 - 0.5e0*v[184] * v[414];
	v[188] = v[184] * (-alphaB[2] + v[190]);
	v[189] = v[184] * (alphaB[1] + v[195]);
	v[191] = v[184] * (alphaB[2] + v[190]);
	v[193] = 1e0 - 0.5e0*v[184] * (v[186] + v[192]);
	v[194] = v[184] * (-alphaB[0] + v[197]);
	v[196] = v[184] * (-alphaB[1] + v[195]);
	v[198] = v[184] * (alphaB[0] + v[197]);
	v[199] = 1e0 - 0.5e0*v[184] * (v[185] + v[192]);
	v[200] = QBi[0][0] * v[187] + QBi[1][0] * v[188] + QBi[2][0] * v[189];
	v[201] = QBi[0][1] * v[187] + QBi[1][1] * v[188] + QBi[2][1] * v[189];
	v[202] = QBi[0][2] * v[187] + QBi[1][2] * v[188] + QBi[2][2] * v[189];
	v[224] = dGBp[0][1] * v[200] + dGBp[1][1] * v[201] + dGBp[2][1] * v[202];
	v[221] = dGBp[0][0] * v[200] + dGBp[1][0] * v[201] + dGBp[2][0] * v[202];
	v[203] = QBi[0][0] * v[191] + QBi[1][0] * v[193] + QBi[2][0] * v[194];
	v[204] = QBi[0][1] * v[191] + QBi[1][1] * v[193] + QBi[2][1] * v[194];
	v[205] = QBi[0][2] * v[191] + QBi[1][2] * v[193] + QBi[2][2] * v[194];
	v[225] = dGBp[0][1] * v[203] + dGBp[1][1] * v[204] + dGBp[2][1] * v[205];
	v[222] = dGBp[0][0] * v[203] + dGBp[1][0] * v[204] + dGBp[2][0] * v[205];
	v[241] = -(v[222] * v[224]) + v[221] * v[225];
	v[206] = QBi[0][0] * v[196] + QBi[1][0] * v[198] + QBi[2][0] * v[199];
	v[207] = QBi[0][1] * v[196] + QBi[1][1] * v[198] + QBi[2][1] * v[199];
	v[208] = QBi[0][2] * v[196] + QBi[1][2] * v[198] + QBi[2][2] * v[199];
	v[226] = dGBp[0][1] * v[206] + dGBp[1][1] * v[207] + dGBp[2][1] * v[208];
	v[223] = dGBp[0][0] * v[206] + dGBp[1][0] * v[207] + dGBp[2][0] * v[208];
	v[239] = v[223] * v[224] - v[221] * v[226];
	v[248] = uA[0] - uB[0] + GAp[0] * v[169] + GAp[1] * v[170] + GAp[2] * v[171] - GBp[0] * v[200] - GBp[1] * v[201]
		- GBp[2] * v[202] + xAi[0] - xBi[0];
	v[249] = uA[1] - uB[1] + GAp[0] * v[172] + GAp[1] * v[173] + GAp[2] * v[174] - GBp[0] * v[203] - GBp[1] * v[204]
		- GBp[2] * v[205] + xAi[1] - xBi[1];
	v[250] = uA[2] - uB[2] + GAp[0] * v[175] + GAp[1] * v[176] + GAp[2] * v[177] - GBp[0] * v[206] - GBp[1] * v[207]
		- GBp[2] * v[208] + xAi[2] - xBi[2];
	v[227] = ((*invertnormalA) ? -1 : 1);
	v[228] = ((*invertnormalB) ? -1 : 1);
	v[229] = -(v[217] * v[219]) + v[216] * v[220];
	v[270] = (v[229] * v[229]) + (v[232] * v[232]) + (v[234] * v[234]);
	v[451] = 1e0 / Power(v[270], 2);
	v[231] = 1e0 / sqrt(v[270]);
	v[236] = -(v[223] * v[225]) + v[222] * v[226];
	v[268] = (v[236] * v[236]) + (v[239] * v[239]) + (v[241] * v[241]);
	v[450] = 1e0 / Power(v[268], 2);
	v[238] = 1e0 / sqrt(v[268]);
	if ((*fixnormal)) {
		v[419] = -normalA[0] + normalB[0];
		v[418] = -normalA[1] + normalB[1];
		v[417] = -normalA[2] + normalB[2];
		v[247] = (v[250] * v[417] + v[249] * v[418] + v[248] * v[419]) / 2e0;
	}
	else {
		v[416] = v[227] * v[231];
		v[415] = v[228] * v[238];
		v[420] = v[241] * v[415] - v[234] * v[416];
		v[421] = v[239] * v[415] - v[232] * v[416];
		v[422] = v[236] * v[415] - v[229] * v[416];
		v[247] = (v[250] * v[420] + v[249] * v[421] + v[248] * v[422]) / 2e0;
	};
	b254 = (*fixnormal);
	if (b254) {
		v[255] = 0e0;
		v[256] = 0e0;
		v[257] = 0e0;
		v[258] = 0e0;
		v[259] = 0e0;
		v[260] = 0e0;
		v[261] = 0e0;
		v[262] = 0e0;
		v[263] = v[417] / 2e0;
		v[264] = v[418] / 2e0;
		v[265] = v[419] / 2e0;
	}
	else {
		v[267] = -v[416] / 2e0;
		v[266] = v[415] / 2e0;
		v[260] = (v[228] * (v[236] * v[248] + v[239] * v[249] + v[241] * v[250])) / 2e0;
		v[259] = v[248] * v[266];
		v[256] = -(v[227] * (v[229] * v[248] + v[232] * v[249] + v[234] * v[250])) / 2e0;
		v[255] = v[248] * v[267];
		v[263] = v[420] / 2e0;
		v[264] = v[421] / 2e0;
		v[265] = v[422] / 2e0;
		v[261] = v[249] * v[266];
		v[262] = v[250] * v[266];
		v[257] = v[249] * v[267];
		v[258] = v[250] * v[267];
	};
	v[424] = -(v[256] / v[270]);
	v[271] = v[231] * v[424];
	v[423] = -(v[260] / v[268]);
	v[269] = v[238] * v[423];
	v[259] = v[259] + v[236] * v[269];
	v[261] = v[261] + v[239] * v[269];
	v[262] = v[262] + v[241] * v[269];
	v[255] = v[255] + v[229] * v[271];
	v[257] = v[257] + v[232] * v[271];
	v[258] = v[258] + v[234] * v[271];
	v[272] = -(v[208] * v[263]) - v[205] * v[264] - v[202] * v[265];
	v[273] = -(v[207] * v[263]) - v[204] * v[264] - v[201] * v[265];
	v[274] = -(v[206] * v[263]) - v[203] * v[264] - v[200] * v[265];
	v[275] = v[177] * v[263] + v[174] * v[264] + v[171] * v[265];
	v[276] = v[176] * v[263] + v[173] * v[264] + v[170] * v[265];
	v[277] = v[175] * v[263] + v[172] * v[264] + v[169] * v[265];
	v[278] = -(v[225] * v[259]) + v[224] * v[261];
	v[279] = v[222] * v[259] - v[221] * v[261];
	v[280] = v[226] * v[259] - v[224] * v[262];
	v[281] = -(v[223] * v[259]) + v[221] * v[262];
	v[282] = -(v[226] * v[261]) + v[225] * v[262];
	v[283] = v[223] * v[261] - v[222] * v[262];
	v[284] = v[208] * v[278] + v[205] * v[280] + v[202] * v[282];
	v[285] = v[207] * v[278] + v[204] * v[280] + v[201] * v[282];
	v[286] = v[206] * v[278] + v[203] * v[280] + v[200] * v[282];
	v[287] = v[208] * v[279] + v[205] * v[281] + v[202] * v[283];
	v[288] = v[207] * v[279] + v[204] * v[281] + v[201] * v[283];
	v[289] = v[206] * v[279] + v[203] * v[281] + v[200] * v[283];
	v[636] = 0e0;
	v[637] = 0e0;
	v[638] = dddGBp[2][0][0][1] * v[284] + dddGBp[1][0][0][1] * v[285] + dddGBp[0][0][0][1] * v[286]
		+ dddGBp[2][1][0][1] * v[287] + dddGBp[1][1][0][1] * v[288] + dddGBp[0][1][0][1] * v[289];
	v[639] = dddGBp[2][0][1][1] * v[284] + dddGBp[1][0][1][1] * v[285] + dddGBp[0][0][1][1] * v[286]
		+ dddGBp[2][1][1][1] * v[287] + dddGBp[1][1][1][1] * v[288] + dddGBp[0][1][1][1] * v[289];
	v[640] = 0e0;
	v[641] = 0e0;
	v[642] = dddGBp[2][0][0][0] * v[284] + dddGBp[1][0][0][0] * v[285] + dddGBp[0][0][0][0] * v[286]
		+ dddGBp[2][1][0][0] * v[287] + dddGBp[1][1][0][0] * v[288] + dddGBp[0][1][0][0] * v[289];
	v[643] = dddGBp[2][0][1][0] * v[284] + dddGBp[1][0][1][0] * v[285] + dddGBp[0][0][1][0] * v[286]
		+ dddGBp[2][1][1][0] * v[287] + dddGBp[1][1][1][0] * v[288] + dddGBp[0][1][1][0] * v[289];
	v[290] = -(v[219] * v[255]) + v[218] * v[257];
	v[291] = v[216] * v[255] - v[215] * v[257];
	v[292] = v[220] * v[255] - v[218] * v[258];
	v[293] = -(v[217] * v[255]) + v[215] * v[258];
	v[294] = -(v[220] * v[257]) + v[219] * v[258];
	v[295] = v[217] * v[257] - v[216] * v[258];
	v[296] = v[177] * v[290] + v[174] * v[292] + v[171] * v[294];
	v[297] = v[176] * v[290] + v[173] * v[292] + v[170] * v[294];
	v[298] = v[175] * v[290] + v[172] * v[292] + v[169] * v[294];
	v[299] = v[177] * v[291] + v[174] * v[293] + v[171] * v[295];
	v[300] = v[176] * v[291] + v[173] * v[293] + v[170] * v[295];
	v[301] = v[175] * v[291] + v[172] * v[293] + v[169] * v[295];
	v[668] = dddGAp[2][0][0][1] * v[296] + dddGAp[1][0][0][1] * v[297] + dddGAp[0][0][0][1] * v[298]
		+ dddGAp[2][1][0][1] * v[299] + dddGAp[1][1][0][1] * v[300] + dddGAp[0][1][0][1] * v[301];
	v[669] = dddGAp[2][0][1][1] * v[296] + dddGAp[1][0][1][1] * v[297] + dddGAp[0][0][1][1] * v[298]
		+ dddGAp[2][1][1][1] * v[299] + dddGAp[1][1][1][1] * v[300] + dddGAp[0][1][1][1] * v[301];
	v[670] = 0e0;
	v[671] = 0e0;
	v[672] = dddGAp[2][0][0][0] * v[296] + dddGAp[1][0][0][0] * v[297] + dddGAp[0][0][0][0] * v[298]
		+ dddGAp[2][1][0][0] * v[299] + dddGAp[1][1][0][0] * v[300] + dddGAp[0][1][0][0] * v[301];
	v[673] = dddGAp[2][0][1][0] * v[296] + dddGAp[1][0][1][0] * v[297] + dddGAp[0][0][1][0] * v[298]
		+ dddGAp[2][1][1][0] * v[299] + dddGAp[1][1][1][0] * v[300] + dddGAp[0][1][1][0] * v[301];
	v[674] = 0e0;
	v[675] = 0e0;
	v[532] = dGAp[2][0] * v[275] + dGAp[1][0] * v[276] + dGAp[0][0] * v[277] + ddGAp[2][0][0] * v[296]
		+ ddGAp[1][0][0] * v[297] + ddGAp[0][0][0] * v[298] + ddGAp[2][1][0] * v[299] + ddGAp[1][1][0] * v[300]
		+ ddGAp[0][1][0] * v[301];
	v[533] = dGAp[2][1] * v[275] + dGAp[1][1] * v[276] + dGAp[0][1] * v[277] + ddGAp[2][0][1] * v[296]
		+ ddGAp[1][0][1] * v[297] + ddGAp[0][0][1] * v[298] + ddGAp[2][1][1] * v[299] + ddGAp[1][1][1] * v[300]
		+ ddGAp[0][1][1] * v[301];
	v[534] = dGBp[2][0] * v[272] + dGBp[1][0] * v[273] + dGBp[0][0] * v[274] + ddGBp[2][0][0] * v[284]
		+ ddGBp[1][0][0] * v[285] + ddGBp[0][0][0] * v[286] + ddGBp[2][1][0] * v[287] + ddGBp[1][1][0] * v[288]
		+ ddGBp[0][1][0] * v[289];
	v[535] = dGBp[2][1] * v[272] + dGBp[1][1] * v[273] + dGBp[0][1] * v[274] + ddGBp[2][0][1] * v[284]
		+ ddGBp[1][0][1] * v[285] + ddGBp[0][0][1] * v[286] + ddGBp[2][1][1] * v[287] + ddGBp[1][1][1] * v[288]
		+ ddGBp[0][1][1] * v[289];
	
	for (i252 = 1; i252 <= 4; i252++) {
		i428 = (i252 == 2 ? 1 : 0);
		i427 = (i252 == 1 ? 1 : 0);
		i426 = (i252 == 4 ? 1 : 0);
		i425 = (i252 == 3 ? 1 : 0);
		v[308] = v[539 + i252];
		v[309] = v[543 + i252];
		v[310] = v[547 + i252];
		v[311] = v[551 + i252];
		v[312] = v[555 + i252];
		v[313] = v[559 + i252];
		v[314] = v[563 + i252];
		v[315] = v[567 + i252];
		v[316] = v[571 + i252];
		v[317] = v[575 + i252];
		v[318] = v[579 + i252];
		v[319] = v[583 + i252];
		v[326] = v[169] * v[308] + v[170] * v[309] + v[171] * v[310];
		v[327] = v[172] * v[308] + v[173] * v[309] + v[174] * v[310];
		v[328] = v[175] * v[308] + v[176] * v[309] + v[177] * v[310];
		v[329] = v[169] * v[311] + v[170] * v[312] + v[171] * v[313];
		v[330] = v[172] * v[311] + v[173] * v[312] + v[174] * v[313];
		v[342] = -(v[216] * v[326]) + v[215] * v[327] + v[219] * v[329] - v[218] * v[330];
		v[331] = v[175] * v[311] + v[176] * v[312] + v[177] * v[313];
		v[348] = -(v[217] * v[327]) + v[216] * v[328] + v[220] * v[330] - v[219] * v[331];
		v[345] = v[217] * v[326] - v[215] * v[328] - v[220] * v[329] + v[218] * v[331];
		v[332] = v[200] * v[314] + v[201] * v[315] + v[202] * v[316];
		v[333] = v[203] * v[314] + v[204] * v[315] + v[205] * v[316];
		v[334] = v[206] * v[314] + v[207] * v[315] + v[208] * v[316];
		v[335] = v[200] * v[317] + v[201] * v[318] + v[202] * v[319];
		v[336] = v[203] * v[317] + v[204] * v[318] + v[205] * v[319];
		v[352] = -(v[222] * v[332]) + v[221] * v[333] + v[225] * v[335] - v[224] * v[336];
		v[337] = v[206] * v[317] + v[207] * v[318] + v[208] * v[319];
		v[358] = -(v[223] * v[333]) + v[222] * v[334] + v[226] * v[336] - v[225] * v[337];
		v[355] = v[223] * v[332] - v[221] * v[334] - v[226] * v[335] + v[224] * v[337];
		v[341] = v[271] * v[342];
		v[344] = v[271] * v[345];
		v[347] = v[234] * v[342] + v[232] * v[345] + v[229] * v[348];
		v[349] = v[271] * v[348];
		v[351] = v[269] * v[352];
		v[354] = v[269] * v[355];
		v[357] = v[241] * v[352] + v[239] * v[355] + v[236] * v[358];
		v[359] = v[269] * v[358];
		v[362] = v[357] * v[423];
		v[364] = v[347] * v[424];
		b365 = (*fixnormal);
		if (b365) {
			v[366] = 0e0;
			v[367] = 0e0;
			v[368] = 0e0;
		}
		else {
			v[325] = v[607 + i252];
			v[324] = v[603 + i252];
			v[323] = v[599 + i252];
			v[322] = v[595 + i252];
			v[321] = v[591 + i252];
			v[320] = v[587 + i252];
			v[340] = v[169] * v[320] + v[170] * v[321] + v[171] * v[322] - v[200] * v[323] - v[201] * v[324] - v[202] * v[325];
			v[339] = v[172] * v[320] + v[173] * v[321] + v[174] * v[322] - v[203] * v[323] - v[204] * v[324] - v[205] * v[325];
			v[338] = v[175] * v[320] + v[176] * v[321] + v[177] * v[322] - v[206] * v[323] - v[207] * v[324] - v[208] * v[325];
			v[370] = -(v[357] * v[415]) / (2e0*v[268]);
			v[369] = (v[347] * v[416]) / (2e0*v[270]);
			v[349] = v[267] * v[340] + v[349] + v[248] * v[369];
			v[344] = v[267] * v[339] + v[344] + v[249] * v[369];
			v[341] = v[267] * v[338] + v[341] + v[250] * v[369];
			v[359] = v[266] * v[340] + v[359] + v[248] * v[370];
			v[368] = v[267] * v[342] + v[266] * v[352] + v[234] * v[369] + v[241] * v[370];
			v[367] = v[267] * v[345] + v[266] * v[355] + v[232] * v[369] + v[239] * v[370];
			v[366] = v[267] * v[348] + v[266] * v[358] + v[229] * v[369] + v[236] * v[370];
			v[354] = v[266] * v[339] + v[354] + v[249] * v[370];
			v[351] = v[266] * v[338] + v[351] + v[250] * v[370];
			v[362] = (v[228] * (v[241] * v[338] + v[239] * v[339] + v[236] * v[340] + v[250] * v[352] + v[249] * v[355]
				+ v[248] * v[358])) / 2e0 + v[362];
			v[364] = -(v[227] * (v[234] * v[338] + v[232] * v[339] + v[229] * v[340] + v[250] * v[342] + v[249] * v[345]
				+ v[248] * v[348])) / 2e0 + v[364];
		};
		v[371] = -(v[238] * (-2e0*v[260] * v[357] + v[268] * v[362])*v[450]) / 2e0;
		v[359] = v[359] + 2e0*v[236] * v[371];
		v[354] = v[354] + 2e0*v[239] * v[371];
		v[351] = v[351] + 2e0*v[241] * v[371];
		v[372] = -(v[231] * (-2e0*v[256] * v[347] + v[270] * v[364])*v[451]) / 2e0;
		v[349] = v[349] + 2e0*v[229] * v[372];
		v[344] = v[344] + 2e0*v[232] * v[372];
		v[341] = v[341] + 2e0*v[234] * v[372];
		v[373] = -(v[202] * v[366]) - v[205] * v[367] - v[208] * v[368];
		v[374] = -(v[201] * v[366]) - v[204] * v[367] - v[207] * v[368];
		v[375] = -(v[200] * v[366]) - v[203] * v[367] - v[206] * v[368];
		v[376] = v[171] * v[366] + v[174] * v[367] + v[177] * v[368];
		v[377] = v[170] * v[366] + v[173] * v[367] + v[176] * v[368];
		v[378] = v[169] * v[366] + v[172] * v[367] + v[175] * v[368];
		v[379] = v[261] * v[332] - v[259] * v[333] + v[224] * v[354] - v[225] * v[359];
		v[380] = -(v[261] * v[335]) + v[259] * v[336] - v[221] * v[354] + v[222] * v[359];
		v[381] = -(v[262] * v[332]) + v[259] * v[334] - v[224] * v[351] + v[226] * v[359];
		v[382] = v[262] * v[335] - v[259] * v[337] + v[221] * v[351] - v[223] * v[359];
		v[383] = v[262] * v[333] - v[261] * v[334] + v[225] * v[351] - v[226] * v[354];
		v[384] = -(v[262] * v[336]) + v[261] * v[337] - v[222] * v[351] + v[223] * v[354];
		v[385] = i425 * v[272] + v[208] * v[379] + v[205] * v[381] + v[202] * v[383];
		v[387] = i425 * v[273] + v[207] * v[379] + v[204] * v[381] + v[201] * v[383];
		v[388] = i425 * v[274] + v[206] * v[379] + v[203] * v[381] + v[200] * v[383];
		v[389] = i426 * v[272] + v[208] * v[380] + v[205] * v[382] + v[202] * v[384];
		v[391] = i426 * v[273] + v[207] * v[380] + v[204] * v[382] + v[201] * v[384];
		v[392] = i426 * v[274] + v[206] * v[380] + v[203] * v[382] + v[200] * v[384];
		v[393] = v[257] * v[326] - v[255] * v[327] + v[218] * v[344] - v[219] * v[349];
		v[394] = -(v[257] * v[329]) + v[255] * v[330] - v[215] * v[344] + v[216] * v[349];
		v[395] = -(v[258] * v[326]) + v[255] * v[328] - v[218] * v[341] + v[220] * v[349];
		v[396] = v[258] * v[329] - v[255] * v[331] + v[215] * v[341] - v[217] * v[349];
		v[397] = v[258] * v[327] - v[257] * v[328] + v[219] * v[341] - v[220] * v[344];
		v[398] = -(v[258] * v[330]) + v[257] * v[331] - v[216] * v[341] + v[217] * v[344];
		v[399] = i427 * v[275] + v[177] * v[393] + v[174] * v[395] + v[171] * v[397];
		v[401] = i427 * v[276] + v[176] * v[393] + v[173] * v[395] + v[170] * v[397];
		v[402] = i427 * v[277] + v[175] * v[393] + v[172] * v[395] + v[169] * v[397];
		v[403] = i428 * v[275] + v[177] * v[394] + v[174] * v[396] + v[171] * v[398];
		v[405] = i428 * v[276] + v[176] * v[394] + v[173] * v[396] + v[170] * v[398];
		v[406] = i428 * v[277] + v[175] * v[394] + v[172] * v[396] + v[169] * v[398];
		v[676] = dGAp[2][0] * v[376] + dGAp[1][0] * v[377] + dGAp[0][0] * v[378] + ddGAp[2][0][0] * v[399]
			+ ddGAp[1][0][0] * v[401] + ddGAp[0][0][0] * v[402] + ddGAp[2][1][0] * v[403] + ddGAp[1][1][0] * v[405]
			+ ddGAp[0][1][0] * v[406] + v[671 + i252];
		v[677] = dGAp[2][1] * v[376] + dGAp[1][1] * v[377] + dGAp[0][1] * v[378] + ddGAp[2][0][1] * v[399]
			+ ddGAp[1][0][1] * v[401] + ddGAp[0][0][1] * v[402] + ddGAp[2][1][1] * v[403] + ddGAp[1][1][1] * v[405]
			+ ddGAp[0][1][1] * v[406] + v[667 + i252];
		v[678] = dGBp[2][0] * v[373] + dGBp[1][0] * v[374] + dGBp[0][0] * v[375] + ddGBp[2][0][0] * v[385]
			+ ddGBp[1][0][0] * v[387] + ddGBp[0][0][0] * v[388] + ddGBp[2][1][0] * v[389] + ddGBp[1][1][0] * v[391]
			+ ddGBp[0][1][0] * v[392] + v[639 + i252];
		v[679] = dGBp[2][1] * v[373] + dGBp[1][1] * v[374] + dGBp[0][1] * v[375] + ddGBp[2][0][1] * v[385]
			+ ddGBp[1][0][1] * v[387] + ddGBp[0][0][1] * v[388] + ddGBp[2][1][1] * v[389] + ddGBp[1][1][1] * v[391]
			+ ddGBp[0][1][1] * v[392] + v[635 + i252];
		
		for (i305 = i252; i305 <= 4; i305++) {
			v[408] = v[675 + i305];
			Hes[i252 - 1][i305 - 1] = v[408];
			if (i252 != i305) {
				Hes[i305 - 1][i252 - 1] = v[408];
			}
			else {
			};
		};/* end for */
	};/* end for */

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mHes(i, j) = Hes[i][j];
}

//Verifica range de coordenadas convectivas
int RigidNURBS_1_RigidNURBS_1::VerifyConvectiveRange(Matrix& mc)
{
	//Retornos:
	//0 - Range fisico da superficie
	//4 - Fora do range fisico da superficie - proximo
	//2 - Fora do range fisico da superficie - distante (nao strong) 

	NURBSSurface* surf1;		//Ponteiro para a superfície 1
	NURBSSurface* surf2;		//Ponteiro para a superfície 2
	surf1 = static_cast<RigidNURBS_1*>(db.surfaces[surf1_ID - 1])->surf;
	surf2 = static_cast<RigidNURBS_1*>(db.surfaces[surf2_ID - 1])->surf;
	
	//Variaveis de medicao de range de coordenadas convectivas
	double range_u1;
	double range_v1;
	double range_u2;
	double range_v2;

	int return_vector[4];
	for (int i = 0; i < 4; i++)
		return_vector[i] = 0;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	range_u1 = perc * convective_range(0, 0);
	if (mc(0, 0) >= surf1->U_knot_vector[surf1->U_order] && mc(0, 0) <= surf1->U_knot_vector[surf1->U_order + surf1->U_dim])
		return_vector[0] = 0;
	else
	{
		if (mc(0, 0) >= (surf1->U_knot_vector[surf1->U_order]-range_u1) && mc(0, 0) <= (surf1->U_knot_vector[surf1->U_order + surf1->U_dim] + range_u1))
			return_vector[0] = 4;
		else
			return_vector[0] = 2;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	range_v1 = perc * convective_range(1, 0);
	if (mc(1, 0) >= surf1->V_knot_vector[surf1->V_order] && mc(1, 0) <= surf1->V_knot_vector[surf1->V_order + surf1->V_dim])
		return_vector[1] = 0;
	else
	{
		if (mc(1, 0) >= (surf1->V_knot_vector[surf1->V_order] - range_v1) && mc(1, 0) <= (surf1->V_knot_vector[surf1->V_order + surf1->V_dim] + range_v1))
			return_vector[1] = 4;
		else
			return_vector[1] = 2;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	range_u2 = perc * convective_range(2, 0);
	if (mc(2, 0) >= surf2->U_knot_vector[surf2->U_order] && mc(2, 0) <= surf2->U_knot_vector[surf2->U_order + surf2->U_dim])
		return_vector[2] = 0;
	else
	{
		if (mc(2, 0) >= (surf2->U_knot_vector[surf2->U_order] - range_u2) && mc(2, 0) <= (surf2->U_knot_vector[surf2->U_order + surf2->U_dim] + range_u2))
			return_vector[2] = 4;
		else
			return_vector[2] = 2;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	range_v2 = perc * convective_range(3, 0);
	if (mc(3, 0) >= surf2->V_knot_vector[surf2->V_order] && mc(3, 0) <= surf2->V_knot_vector[surf2->V_order + surf2->V_dim])
		return_vector[3] = 0;
	else
	{
		if (mc(3, 0) >= (surf2->V_knot_vector[surf2->V_order] - range_v2) && mc(3, 0) <= (surf2->V_knot_vector[surf2->V_order + surf2->V_dim] + range_v2))
			return_vector[3] = 4;
		else
			return_vector[3] = 2;
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

//Initialize range of validity of convective coordinates
void RigidNURBS_1_RigidNURBS_1::InitializeConvectiveRange()
{
	NURBSSurface* surf1;		//Ponteiro para a superfície 1
	NURBSSurface* surf2;		//Ponteiro para a superfície 2
	surf1 = static_cast<RigidNURBS_1*>(db.surfaces[surf1_ID - 1])->surf;
	surf2 = static_cast<RigidNURBS_1*>(db.surfaces[surf2_ID - 1])->surf;

	RigidNURBS_1* surf1a;		//Ponteiro para a superfície 1
	RigidNURBS_1* surf2a;		//Ponteiro para a superfície 2
	surf1a = static_cast<RigidNURBS_1*>(db.surfaces[surf1_ID - 1]);
	surf2a = static_cast<RigidNURBS_1*>(db.surfaces[surf2_ID - 1]);

	convective_min(0, 0) = surf1->U_knot_vector[surf1->U_order];
	convective_max(0, 0) = surf1->U_knot_vector[surf1->U_order + surf1->U_dim];
	//Degeneration
	if (surf1a->degeneration[0] == true)
	{
		convective_range(0, 0) = 0.0;
		
	}
	//No degeneration
	else
	{
		convective_range(0, 0) = surf1->U_knot_vector[surf1->U_order + surf1->U_dim] - surf1->U_knot_vector[surf1->U_order];
		
	}

	convective_min(1, 0) = surf1->V_knot_vector[surf1->V_order];
	convective_max(1, 0) = surf1->V_knot_vector[surf1->V_order + surf1->V_dim];
	//Degeneration
	if (surf1a->degeneration[1] == true)
	{
		convective_range(1, 0) = 0.0;
		
	}
	//No degeneration
	else
	{
		convective_range(1, 0) = surf1->V_knot_vector[surf1->V_order + surf1->V_dim] - surf1->V_knot_vector[surf1->V_order];
		
	}

	convective_min(2, 0) = surf2->U_knot_vector[surf2->U_order];
	convective_max(2, 0) = surf2->U_knot_vector[surf2->U_order + surf2->U_dim];
	//Degeneration
	if (surf2a->degeneration[0] == true)
	{
		convective_range(2, 0) = 0.0;
		
	}
	//No degeneration
	else
	{
		convective_range(2, 0) = surf2->U_knot_vector[surf2->U_order + surf2->U_dim] - surf2->U_knot_vector[surf2->U_order];
		
	}

	convective_min(3, 0) = surf2->V_knot_vector[surf2->V_order];
	convective_max(3, 0) = surf2->V_knot_vector[surf2->V_order + surf2->V_dim];
	//Degeneration
	if (surf2a->degeneration[1] == true)
	{
		convective_range(3, 0) = 0.0;
		
	}
	//No degeneration
	else
	{
		convective_range(3, 0) = surf2->V_knot_vector[surf2->V_order + surf2->V_dim] - surf2->V_knot_vector[surf2->V_order];
		
	}
}
