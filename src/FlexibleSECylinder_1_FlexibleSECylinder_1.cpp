#include "FlexibleSECylinder_1_FlexibleSECylinder_1.h"

#include"Database.h"
//Variáveis globais
extern
Database db; 

FlexibleSECylinder_1_FlexibleSECylinder_1::FlexibleSECylinder_1_FlexibleSECylinder_1()
{
	DefaultValues();
}

FlexibleSECylinder_1_FlexibleSECylinder_1::~FlexibleSECylinder_1_FlexibleSECylinder_1()
{
	Free();
}

void FlexibleSECylinder_1_FlexibleSECylinder_1::InitializeConvectiveRange()
{
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);

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
		convective_range(0, 0) = 2;
	}

	convective_min(1, 0) = 0;
	convective_max(1, 0) = 2 * PI;
	//Degeneration
	if (surf1->degeneration[1] == true)
	{
		convective_range(1, 0) = 0.0;
	}
	//No degeneration
	else
	{
		convective_range(1, 0) = 2 * PI;
		
	}

	convective_min(2, 0) = -1;
	convective_max(2, 0) = +1;
	//Degeneration
	if (surf2->degeneration[0] == true)
	{
		convective_range(2, 0) = 0.0;
		
	}
	//No degeneration
	else
	{
		convective_range(2, 0) = 2;
		
	}

	convective_min(3, 0) = 0;
	convective_max(3, 0) = 2 * PI;
	//Degeneration
	if (surf2->degeneration[1] == true)
	{
		convective_range(3, 0) = 0.0;
		
	}
	//No degeneration
	else
	{
		convective_range(3, 0) = 2 * PI;
		
	}
}

//Verifica range de coordenadas convectivas
int FlexibleSECylinder_1_FlexibleSECylinder_1::VerifyConvectiveRange(Matrix& mc)
{
	int return_value;
	//Retornos:
	//0 - Range fisico da superficie
	//4 - Fora do range fisico da superficie - proximo
	//2 - Fora do range fisico da superficie - distante

	//Se está no range local de interesse - domínio físico da superfície
	if (abs(mc(0, 0)) <= 1.0 && abs(mc(2, 0)) <= 1.0)
		return_value = 0;	//Houve convergência, está no range físico - forte candidato a contato
	else
	{
		if (abs(mc(0, 0)) <= (1.00+perc*2) && abs(mc(2, 0)) <= (1.0+perc*2))
			return_value = 4;	//Houve convergência, está no range físico - forte candidato a contato
		else
			return_value = 2;	//Houve convergência, mas não está no range físico - deve ser monitorado com cuidado - possivelmente superfície vizinha
	}
	return return_value;
}

//Calcula e rotorna o gap (com sinal)
double FlexibleSECylinder_1_FlexibleSECylinder_1::Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variável temporária - AceGen
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
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
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	aA = &surf1->a;
	aB = &surf2->a;
	bA = &surf1->b;
	bB = &surf2->b;
	eA = &surf1->e;
	eB = &surf2->e;
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_A->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_A->getMatrix();
	dduiB = surf2->ddui_A->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_AAi->getMatrix();
	xBBi = surf2->x_BAi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	QBBi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
		QBBi[i] = new double[3];
	}
	//Salvando variáveis locais para montagem de superfícies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_AAi->MatrixToPtr(QABi, 3);
	surf2->Q_BAi->MatrixToPtr(QBBi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double gap;
	double *c = mc.getMatrix();
	
	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	//AceGen
	int b303;
	v[247] = cos(c[1]);
	v[169] = sin(c[1]);
	v[251] = fabs(v[169]);
	v[249] = fabs(v[247]);
	v[267] = cos(c[3]);
	v[242] = sin(c[3]);
	v[271] = fabs(v[242]);
	v[269] = fabs(v[267]);
	v[109] = Power(dA[3], 2);
	v[107] = (dA[3] * dA[4]) / 2e0;
	v[102] = Power(dA[4], 2);
	v[114] = (dA[4] * dA[5]) / 2e0;
	v[112] = (dA[3] * dA[5]) / 2e0;
	v[103] = Power(dA[5], 2);
	v[312] = v[102] + v[103];
	v[128] = Power(dA[9], 2);
	v[126] = (dA[10] * dA[9]) / 2e0;
	v[121] = Power(dA[10], 2);
	v[133] = (dA[10] * dA[11]) / 2e0;
	v[131] = (dA[11] * dA[9]) / 2e0;
	v[122] = Power(dA[11], 2);
	v[313] = v[121] + v[122];
	v[182] = Power(dB[3], 2);
	v[180] = (dB[3] * dB[4]) / 2e0;
	v[175] = Power(dB[4], 2);
	v[187] = (dB[4] * dB[5]) / 2e0;
	v[185] = (dB[3] * dB[5]) / 2e0;
	v[176] = Power(dB[5], 2);
	v[317] = v[175] + v[176];
	v[201] = Power(dB[9], 2);
	v[199] = (dB[10] * dB[9]) / 2e0;
	v[194] = Power(dB[10], 2);
	v[206] = (dB[10] * dB[11]) / 2e0;
	v[204] = (dB[11] * dB[9]) / 2e0;
	v[195] = Power(dB[11], 2);
	v[318] = v[194] + v[195];
	v[101] = 4e0 / (4e0 + v[109] + v[312]);
	v[104] = 1e0 - (v[101] * v[312]) / 2e0;
	v[105] = v[101] * (-dA[5] + v[107]);
	v[106] = v[101] * (dA[4] + v[112]);
	v[108] = v[101] * (dA[5] + v[107]);
	v[110] = 1e0 - (v[101] * (v[103] + v[109])) / 2e0;
	v[111] = v[101] * (-dA[3] + v[114]);
	v[113] = v[101] * (-dA[4] + v[112]);
	v[115] = v[101] * (dA[3] + v[114]);
	v[116] = 1e0 - (v[101] * (v[102] + v[109])) / 2e0;
	v[120] = 4e0 / (4e0 + v[128] + v[313]);
	v[123] = 1e0 - (v[120] * v[313]) / 2e0;
	v[124] = v[120] * (-dA[11] + v[126]);
	v[125] = v[120] * (dA[10] + v[131]);
	v[127] = v[120] * (dA[11] + v[126]);
	v[129] = 1e0 - (v[120] * (v[122] + v[128])) / 2e0;
	v[130] = v[120] * (-dA[9] + v[133]);
	v[132] = v[120] * (-dA[10] + v[131]);
	v[134] = v[120] * (dA[9] + v[133]);
	v[135] = 1e0 - (v[120] * (v[121] + v[128])) / 2e0;
	v[139] = QAAi[0][0] * v[104] + QAAi[1][0] * v[105] + QAAi[2][0] * v[106];
	v[140] = QAAi[0][1] * v[104] + QAAi[1][1] * v[105] + QAAi[2][1] * v[106];
	v[142] = QAAi[0][0] * v[108] + QAAi[1][0] * v[110] + QAAi[2][0] * v[111];
	v[143] = QAAi[0][1] * v[108] + QAAi[1][1] * v[110] + QAAi[2][1] * v[111];
	v[145] = QAAi[0][0] * v[113] + QAAi[1][0] * v[115] + QAAi[2][0] * v[116];
	v[146] = QAAi[0][1] * v[113] + QAAi[1][1] * v[115] + QAAi[2][1] * v[116];
	v[148] = QBAi[0][0] * v[123] + QBAi[1][0] * v[124] + QBAi[2][0] * v[125];
	v[149] = QBAi[0][1] * v[123] + QBAi[1][1] * v[124] + QBAi[2][1] * v[125];
	v[151] = QBAi[0][0] * v[127] + QBAi[1][0] * v[129] + QBAi[2][0] * v[130];
	v[152] = QBAi[0][1] * v[127] + QBAi[1][1] * v[129] + QBAi[2][1] * v[130];
	v[154] = QBAi[0][0] * v[132] + QBAi[1][0] * v[134] + QBAi[2][0] * v[135];
	v[155] = QBAi[0][1] * v[132] + QBAi[1][1] * v[134] + QBAi[2][1] * v[135];
	v[163] = (1e0 - c[0]) / 2e0;
	v[164] = (1e0 + c[0]) / 2e0;
	v[165] = 2e0 / (*eA);
	v[250] = -1e0 + v[165];
	v[248] = Power(v[249], v[165]) + Power(v[251], v[165]);
	v[252] = Power(v[248], -1e0 - 1e0 / v[165])*(-(v[247] * Power(v[251], v[250])*_copysign(1.e0, v[169]))
		+ v[169] * Power(v[249], v[250])*_copysign(1.e0, v[247]));
	v[166] = 1e0 / Power(v[248], 1 / v[165]);
	v[316] = v[166] * v[169];
	v[315] = v[166] * v[247];
	v[254] = (*bA)*(v[169] * v[252] + v[315]);
	v[253] = (*aA)*(v[247] * v[252] - v[316]);
	v[257] = v[163] * (v[145] * v[253] + v[146] * v[254]) + v[164] * (v[154] * v[253] + v[155] * v[254]);
	v[256] = v[163] * (v[142] * v[253] + v[143] * v[254]) + v[164] * (v[151] * v[253] + v[152] * v[254]);
	v[255] = v[163] * (v[139] * v[253] + v[140] * v[254]) + v[164] * (v[148] * v[253] + v[149] * v[254]);
	v[168] = (*aA)*v[315];
	v[170] = (*bA)*v[316];
	v[265] = dA[8] + v[154] * v[168] + v[155] * v[170] + xBAi[2];
	v[264] = dA[2] + v[145] * v[168] + v[146] * v[170] + xAAi[2];
	v[266] = (-v[264] + v[265]) / 2e0;
	v[262] = dA[7] + v[151] * v[168] + v[152] * v[170] + xBAi[1];
	v[261] = dA[1] + v[142] * v[168] + v[143] * v[170] + xAAi[1];
	v[263] = (-v[261] + v[262]) / 2e0;
	v[259] = dA[6] + v[148] * v[168] + v[149] * v[170] + xBAi[0];
	v[258] = dA[0] + v[139] * v[168] + v[140] * v[170] + xAAi[0];
	v[260] = (-v[258] + v[259]) / 2e0;
	v[174] = 4e0 / (4e0 + v[182] + v[317]);
	v[177] = 1e0 - (v[174] * v[317]) / 2e0;
	v[178] = v[174] * (-dB[5] + v[180]);
	v[179] = v[174] * (dB[4] + v[185]);
	v[181] = v[174] * (dB[5] + v[180]);
	v[183] = 1e0 - (v[174] * (v[176] + v[182])) / 2e0;
	v[184] = v[174] * (-dB[3] + v[187]);
	v[186] = v[174] * (-dB[4] + v[185]);
	v[188] = v[174] * (dB[3] + v[187]);
	v[189] = 1e0 - (v[174] * (v[175] + v[182])) / 2e0;
	v[193] = 4e0 / (4e0 + v[201] + v[318]);
	v[196] = 1e0 - (v[193] * v[318]) / 2e0;
	v[197] = v[193] * (-dB[11] + v[199]);
	v[198] = v[193] * (dB[10] + v[204]);
	v[200] = v[193] * (dB[11] + v[199]);
	v[202] = 1e0 - (v[193] * (v[195] + v[201])) / 2e0;
	v[203] = v[193] * (-dB[9] + v[206]);
	v[205] = v[193] * (-dB[10] + v[204]);
	v[207] = v[193] * (dB[9] + v[206]);
	v[208] = 1e0 - (v[193] * (v[194] + v[201])) / 2e0;
	v[212] = QABi[0][0] * v[177] + QABi[1][0] * v[178] + QABi[2][0] * v[179];
	v[213] = QABi[0][1] * v[177] + QABi[1][1] * v[178] + QABi[2][1] * v[179];
	v[215] = QABi[0][0] * v[181] + QABi[1][0] * v[183] + QABi[2][0] * v[184];
	v[216] = QABi[0][1] * v[181] + QABi[1][1] * v[183] + QABi[2][1] * v[184];
	v[218] = QABi[0][0] * v[186] + QABi[1][0] * v[188] + QABi[2][0] * v[189];
	v[219] = QABi[0][1] * v[186] + QABi[1][1] * v[188] + QABi[2][1] * v[189];
	v[221] = QBBi[0][0] * v[196] + QBBi[1][0] * v[197] + QBBi[2][0] * v[198];
	v[222] = QBBi[0][1] * v[196] + QBBi[1][1] * v[197] + QBBi[2][1] * v[198];
	v[224] = QBBi[0][0] * v[200] + QBBi[1][0] * v[202] + QBBi[2][0] * v[203];
	v[225] = QBBi[0][1] * v[200] + QBBi[1][1] * v[202] + QBBi[2][1] * v[203];
	v[227] = QBBi[0][0] * v[205] + QBBi[1][0] * v[207] + QBBi[2][0] * v[208];
	v[228] = QBBi[0][1] * v[205] + QBBi[1][1] * v[207] + QBBi[2][1] * v[208];
	v[236] = (1e0 - c[2]) / 2e0;
	v[237] = (1e0 + c[2]) / 2e0;
	v[238] = 2e0 / (*eB);
	v[270] = -1e0 + v[238];
	v[268] = Power(v[269], v[238]) + Power(v[271], v[238]);
	v[272] = Power(v[268], -1e0 - 1e0 / v[238])*(-(v[267] * Power(v[271], v[270])*_copysign(1.e0, v[242]))
		+ v[242] * Power(v[269], v[270])*_copysign(1.e0, v[267]));
	v[239] = 1e0 / Power(v[268], 1 / v[238]);
	v[321] = v[239] * v[242];
	v[320] = v[239] * v[267];
	v[274] = (*bB)*(v[242] * v[272] + v[320]);
	v[273] = (*aB)*(v[267] * v[272] - v[321]);
	v[277] = v[236] * (v[218] * v[273] + v[219] * v[274]) + v[237] * (v[227] * v[273] + v[228] * v[274]);
	v[276] = v[236] * (v[215] * v[273] + v[216] * v[274]) + v[237] * (v[224] * v[273] + v[225] * v[274]);
	v[275] = v[236] * (v[212] * v[273] + v[213] * v[274]) + v[237] * (v[221] * v[273] + v[222] * v[274]);
	v[241] = (*aB)*v[320];
	v[243] = (*bB)*v[321];
	v[285] = dB[8] + v[227] * v[241] + v[228] * v[243] + xBBi[2];
	v[284] = dB[2] + v[218] * v[241] + v[219] * v[243] + xABi[2];
	v[286] = (-v[284] + v[285]) / 2e0;
	v[282] = dB[7] + v[224] * v[241] + v[225] * v[243] + xBBi[1];
	v[281] = dB[1] + v[215] * v[241] + v[216] * v[243] + xABi[1];
	v[283] = (-v[281] + v[282]) / 2e0;
	v[279] = dB[6] + v[221] * v[241] + v[222] * v[243] + xBBi[0];
	v[278] = dB[0] + v[212] * v[241] + v[213] * v[243] + xABi[0];
	v[280] = (-v[278] + v[279]) / 2e0;
	v[308] = v[163] * v[258] + v[164] * v[259] - v[236] * v[278] - v[237] * v[279];
	v[309] = v[163] * v[261] + v[164] * v[262] - v[236] * v[281] - v[237] * v[282];
	v[310] = v[163] * v[264] + v[164] * v[265] - v[236] * v[284] - v[237] * v[285];
	v[289] = -(v[257] * v[263]) + v[256] * v[266];
	v[296] = -(v[277] * v[283]) + v[276] * v[286];
	if ((*fixnormal)) {
		v[307] = 0.5e0*((-normalA[0] + normalB[0])*v[308] + (-normalA[1] + normalB[1])*v[309] + (-normalA[2]
			+ normalB[2])*v[310]);
	}
	else {
		v[301] = -(v[276] * v[280]) + v[275] * v[283];
		v[299] = v[277] * v[280] - v[275] * v[286];
		v[294] = -(v[256] * v[260]) + v[255] * v[263];
		v[292] = v[257] * v[260] - v[255] * v[266];
		v[307] = (-0.5e0*((*normalintA) ? -1 : 1)*(v[289] * v[308] + v[292] * v[309] + v[294] * v[310])) / sqrt(
			(v[289] * v[289]) + (v[292] * v[292]) + (v[294] * v[294])) + (0.5e0*((*normalintB) ? -1 : 1)*(v[296] * v[308]
				+ v[299] * v[309] + v[301] * v[310])) / sqrt((v[296] * v[296]) + (v[299] * v[299]) + (v[301] * v[301]));
	};
	(gap) = v[307];
	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
		delete[] QBBi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
	delete[] QBBi;

	return gap;
}

//Calcula o Gradiente do gap
void FlexibleSECylinder_1_FlexibleSECylinder_1::GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variável temporária - AceGen
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
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
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	aA = &surf1->a;
	aB = &surf2->a;
	bA = &surf1->b;
	bB = &surf2->b;
	eA = &surf1->e;
	eB = &surf2->e;
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_A->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_A->getMatrix();
	dduiB = surf2->ddui_A->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_AAi->getMatrix();
	xBBi = surf2->x_BAi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	QBBi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
		QBBi[i] = new double[3];
	}
	//Salvando variáveis locais para montagem de superfícies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_AAi->MatrixToPtr(QABi, 3);
	surf2->Q_BAi->MatrixToPtr(QBBi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double *c = mc.getMatrix();
	double Gra[4];
	//Pontos de singularidade - função phi
	double test1 = 2 * c[1] / PI;
	test1 = abs(test1 - int(test1));//testa o quão próximo de inteiro é esse número
	if (test1 < tol_small_1)
		c[1] += tol_small_1;
	double test3 = 2 * c[3] / PI;
	test3 = abs(test3 - int(test3));//testa o quão próximo de inteiro é esse número
	if (test3 < tol_small_1)
		c[3] += tol_small_1;

	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	//AceGen
	int i312, i426, i427, b303, b314;
	v[247] = cos(c[1]);
	v[169] = sin(c[1]);
	v[388] = _copysign(1.e0, v[169]);
	v[251] = fabs(v[169]);
	v[391] = _copysign(1.e0, v[247]);
	v[249] = fabs(v[247]);
	v[267] = cos(c[3]);
	v[242] = sin(c[3]);
	v[401] = _copysign(1.e0, v[242]);
	v[271] = fabs(v[242]);
	v[404] = _copysign(1.e0, v[267]);
	v[269] = fabs(v[267]);
	v[109] = Power(dA[3], 2);
	v[107] = (dA[3] * dA[4]) / 2e0;
	v[102] = Power(dA[4], 2);
	v[114] = (dA[4] * dA[5]) / 2e0;
	v[112] = (dA[3] * dA[5]) / 2e0;
	v[103] = Power(dA[5], 2);
	v[411] = v[102] + v[103];
	v[128] = Power(dA[9], 2);
	v[126] = (dA[10] * dA[9]) / 2e0;
	v[121] = Power(dA[10], 2);
	v[133] = (dA[10] * dA[11]) / 2e0;
	v[131] = (dA[11] * dA[9]) / 2e0;
	v[122] = Power(dA[11], 2);
	v[412] = v[121] + v[122];
	v[182] = Power(dB[3], 2);
	v[180] = (dB[3] * dB[4]) / 2e0;
	v[175] = Power(dB[4], 2);
	v[187] = (dB[4] * dB[5]) / 2e0;
	v[185] = (dB[3] * dB[5]) / 2e0;
	v[176] = Power(dB[5], 2);
	v[416] = v[175] + v[176];
	v[201] = Power(dB[9], 2);
	v[199] = (dB[10] * dB[9]) / 2e0;
	v[194] = Power(dB[10], 2);
	v[206] = (dB[10] * dB[11]) / 2e0;
	v[204] = (dB[11] * dB[9]) / 2e0;
	v[195] = Power(dB[11], 2);
	v[417] = v[194] + v[195];
	v[101] = 4e0 / (4e0 + v[109] + v[411]);
	v[104] = 1e0 - (v[101] * v[411]) / 2e0;
	v[105] = v[101] * (-dA[5] + v[107]);
	v[106] = v[101] * (dA[4] + v[112]);
	v[108] = v[101] * (dA[5] + v[107]);
	v[110] = 1e0 - (v[101] * (v[103] + v[109])) / 2e0;
	v[111] = v[101] * (-dA[3] + v[114]);
	v[113] = v[101] * (-dA[4] + v[112]);
	v[115] = v[101] * (dA[3] + v[114]);
	v[116] = 1e0 - (v[101] * (v[102] + v[109])) / 2e0;
	v[120] = 4e0 / (4e0 + v[128] + v[412]);
	v[123] = 1e0 - (v[120] * v[412]) / 2e0;
	v[124] = v[120] * (-dA[11] + v[126]);
	v[125] = v[120] * (dA[10] + v[131]);
	v[127] = v[120] * (dA[11] + v[126]);
	v[129] = 1e0 - (v[120] * (v[122] + v[128])) / 2e0;
	v[130] = v[120] * (-dA[9] + v[133]);
	v[132] = v[120] * (-dA[10] + v[131]);
	v[134] = v[120] * (dA[9] + v[133]);
	v[135] = 1e0 - (v[120] * (v[121] + v[128])) / 2e0;
	v[139] = QAAi[0][0] * v[104] + QAAi[1][0] * v[105] + QAAi[2][0] * v[106];
	v[140] = QAAi[0][1] * v[104] + QAAi[1][1] * v[105] + QAAi[2][1] * v[106];
	v[142] = QAAi[0][0] * v[108] + QAAi[1][0] * v[110] + QAAi[2][0] * v[111];
	v[143] = QAAi[0][1] * v[108] + QAAi[1][1] * v[110] + QAAi[2][1] * v[111];
	v[145] = QAAi[0][0] * v[113] + QAAi[1][0] * v[115] + QAAi[2][0] * v[116];
	v[146] = QAAi[0][1] * v[113] + QAAi[1][1] * v[115] + QAAi[2][1] * v[116];
	v[148] = QBAi[0][0] * v[123] + QBAi[1][0] * v[124] + QBAi[2][0] * v[125];
	v[149] = QBAi[0][1] * v[123] + QBAi[1][1] * v[124] + QBAi[2][1] * v[125];
	v[151] = QBAi[0][0] * v[127] + QBAi[1][0] * v[129] + QBAi[2][0] * v[130];
	v[152] = QBAi[0][1] * v[127] + QBAi[1][1] * v[129] + QBAi[2][1] * v[130];
	v[154] = QBAi[0][0] * v[132] + QBAi[1][0] * v[134] + QBAi[2][0] * v[135];
	v[155] = QBAi[0][1] * v[132] + QBAi[1][1] * v[134] + QBAi[2][1] * v[135];
	v[163] = (1e0 - c[0]) / 2e0;
	v[164] = (1e0 + c[0]) / 2e0;
	v[165] = 2e0 / (*eA);
	v[378] = -1e0 - 1e0 / v[165];
	v[250] = -1e0 + v[165];
	v[394] = -1e0 + v[250];
	v[392] = Power(v[249], v[250]);
	v[390] = Power(v[251], v[250]);
	v[380] = v[165] * (v[247] * v[388] * v[390] - v[169] * v[391] * v[392]);
	v[248] = Power(v[249], v[165]) + Power(v[251], v[165]);
	v[377] = Power(v[248], v[378]);
	v[252] = -((v[377] * v[380]) / v[165]);
	v[166] = 1e0 / Power(v[248], 1 / v[165]);
	v[415] = v[166] * v[169];
	v[414] = v[166] * v[247];
	v[254] = (*bA)*(v[169] * v[252] + v[414]);
	v[253] = (*aA)*(v[247] * v[252] - v[415]);
	v[387] = v[151] * v[253] + v[152] * v[254];
	v[386] = v[148] * v[253] + v[149] * v[254];
	v[385] = v[154] * v[253] + v[155] * v[254];
	v[384] = v[142] * v[253] + v[143] * v[254];
	v[383] = v[139] * v[253] + v[140] * v[254];
	v[382] = v[145] * v[253] + v[146] * v[254];
	v[257] = v[163] * v[382] + v[164] * v[385];
	v[256] = v[163] * v[384] + v[164] * v[387];
	v[255] = v[163] * v[383] + v[164] * v[386];
	v[168] = (*aA)*v[414];
	v[170] = (*bA)*v[415];
	v[265] = dA[8] + v[154] * v[168] + v[155] * v[170] + xBAi[2];
	v[264] = dA[2] + v[145] * v[168] + v[146] * v[170] + xAAi[2];
	v[266] = (-v[264] + v[265]) / 2e0;
	v[262] = dA[7] + v[151] * v[168] + v[152] * v[170] + xBAi[1];
	v[261] = dA[1] + v[142] * v[168] + v[143] * v[170] + xAAi[1];
	v[263] = (-v[261] + v[262]) / 2e0;
	v[259] = dA[6] + v[148] * v[168] + v[149] * v[170] + xBAi[0];
	v[258] = dA[0] + v[139] * v[168] + v[140] * v[170] + xAAi[0];
	v[260] = (-v[258] + v[259]) / 2e0;
	v[294] = -(v[256] * v[260]) + v[255] * v[263];
	v[292] = v[257] * v[260] - v[255] * v[266];
	v[174] = 4e0 / (4e0 + v[182] + v[416]);
	v[177] = 1e0 - (v[174] * v[416]) / 2e0;
	v[178] = v[174] * (-dB[5] + v[180]);
	v[179] = v[174] * (dB[4] + v[185]);
	v[181] = v[174] * (dB[5] + v[180]);
	v[183] = 1e0 - (v[174] * (v[176] + v[182])) / 2e0;
	v[184] = v[174] * (-dB[3] + v[187]);
	v[186] = v[174] * (-dB[4] + v[185]);
	v[188] = v[174] * (dB[3] + v[187]);
	v[189] = 1e0 - (v[174] * (v[175] + v[182])) / 2e0;
	v[193] = 4e0 / (4e0 + v[201] + v[417]);
	v[196] = 1e0 - (v[193] * v[417]) / 2e0;
	v[197] = v[193] * (-dB[11] + v[199]);
	v[198] = v[193] * (dB[10] + v[204]);
	v[200] = v[193] * (dB[11] + v[199]);
	v[202] = 1e0 - (v[193] * (v[195] + v[201])) / 2e0;
	v[203] = v[193] * (-dB[9] + v[206]);
	v[205] = v[193] * (-dB[10] + v[204]);
	v[207] = v[193] * (dB[9] + v[206]);
	v[208] = 1e0 - (v[193] * (v[194] + v[201])) / 2e0;
	v[212] = QABi[0][0] * v[177] + QABi[1][0] * v[178] + QABi[2][0] * v[179];
	v[213] = QABi[0][1] * v[177] + QABi[1][1] * v[178] + QABi[2][1] * v[179];
	v[215] = QABi[0][0] * v[181] + QABi[1][0] * v[183] + QABi[2][0] * v[184];
	v[216] = QABi[0][1] * v[181] + QABi[1][1] * v[183] + QABi[2][1] * v[184];
	v[218] = QABi[0][0] * v[186] + QABi[1][0] * v[188] + QABi[2][0] * v[189];
	v[219] = QABi[0][1] * v[186] + QABi[1][1] * v[188] + QABi[2][1] * v[189];
	v[221] = QBBi[0][0] * v[196] + QBBi[1][0] * v[197] + QBBi[2][0] * v[198];
	v[222] = QBBi[0][1] * v[196] + QBBi[1][1] * v[197] + QBBi[2][1] * v[198];
	v[224] = QBBi[0][0] * v[200] + QBBi[1][0] * v[202] + QBBi[2][0] * v[203];
	v[225] = QBBi[0][1] * v[200] + QBBi[1][1] * v[202] + QBBi[2][1] * v[203];
	v[227] = QBBi[0][0] * v[205] + QBBi[1][0] * v[207] + QBBi[2][0] * v[208];
	v[228] = QBBi[0][1] * v[205] + QBBi[1][1] * v[207] + QBBi[2][1] * v[208];
	v[236] = (1e0 - c[2]) / 2e0;
	v[237] = (1e0 + c[2]) / 2e0;
	v[238] = 2e0 / (*eB);
	v[353] = -1e0 - 1e0 / v[238];
	v[270] = -1e0 + v[238];
	v[407] = -1e0 + v[270];
	v[405] = Power(v[269], v[270]);
	v[403] = Power(v[271], v[270]);
	v[355] = v[238] * (v[267] * v[401] * v[403] - v[242] * v[404] * v[405]);
	v[268] = Power(v[269], v[238]) + Power(v[271], v[238]);
	v[352] = Power(v[268], v[353]);
	v[272] = -((v[352] * v[355]) / v[238]);
	v[239] = 1e0 / Power(v[268], 1 / v[238]);
	v[420] = v[239] * v[242];
	v[419] = v[239] * v[267];
	v[274] = (*bB)*(v[242] * v[272] + v[419]);
	v[273] = (*aB)*(v[267] * v[272] - v[420]);
	v[400] = v[224] * v[273] + v[225] * v[274];
	v[399] = v[221] * v[273] + v[222] * v[274];
	v[398] = v[227] * v[273] + v[228] * v[274];
	v[397] = v[215] * v[273] + v[216] * v[274];
	v[396] = v[212] * v[273] + v[213] * v[274];
	v[395] = v[218] * v[273] + v[219] * v[274];
	v[277] = v[236] * v[395] + v[237] * v[398];
	v[276] = v[236] * v[397] + v[237] * v[400];
	v[275] = v[236] * v[396] + v[237] * v[399];
	v[241] = (*aB)*v[419];
	v[243] = (*bB)*v[420];
	v[285] = dB[8] + v[227] * v[241] + v[228] * v[243] + xBBi[2];
	v[284] = dB[2] + v[218] * v[241] + v[219] * v[243] + xABi[2];
	v[286] = (-v[284] + v[285]) / 2e0;
	v[282] = dB[7] + v[224] * v[241] + v[225] * v[243] + xBBi[1];
	v[281] = dB[1] + v[215] * v[241] + v[216] * v[243] + xABi[1];
	v[283] = (-v[281] + v[282]) / 2e0;
	v[279] = dB[6] + v[221] * v[241] + v[222] * v[243] + xBBi[0];
	v[278] = dB[0] + v[212] * v[241] + v[213] * v[243] + xABi[0];
	v[280] = (-v[278] + v[279]) / 2e0;
	v[301] = -(v[276] * v[280]) + v[275] * v[283];
	v[299] = v[277] * v[280] - v[275] * v[286];
	v[308] = v[163] * v[258] + v[164] * v[259] - v[236] * v[278] - v[237] * v[279];
	v[309] = v[163] * v[261] + v[164] * v[262] - v[236] * v[281] - v[237] * v[282];
	v[310] = v[163] * v[264] + v[164] * v[265] - v[236] * v[284] - v[237] * v[285];
	v[289] = -(v[257] * v[263]) + v[256] * v[266];
	v[330] = (v[289] * v[289]) + (v[292] * v[292]) + (v[294] * v[294]);
	v[291] = 1e0 / sqrt(v[330]);
	v[296] = -(v[277] * v[283]) + v[276] * v[286];
	v[328] = (v[296] * v[296]) + (v[299] * v[299]) + (v[301] * v[301]);
	v[298] = 1e0 / sqrt(v[328]);
	if ((*fixnormal)) {
		v[425] = -normalA[0] + normalB[0];
		v[424] = -normalA[1] + normalB[1];
		v[423] = -normalA[2] + normalB[2];
		v[307] = 0.5e0*(v[310] * v[423] + v[309] * v[424] + v[308] * v[425]);
	}
	else {
		i427 = ((*normalintA) ? -1 : 1);
		i426 = ((*normalintB) ? -1 : 1);
		v[422] = i427 * v[291];
		v[421] = i426 * v[298];
		v[428] = v[301] * v[421] - v[294] * v[422];
		v[429] = v[299] * v[421] - v[292] * v[422];
		v[430] = v[296] * v[421] - v[289] * v[422];
		v[307] = 0.5e0*(v[310] * v[428] + v[309] * v[429] + v[308] * v[430]);
	};
	b314 = (*fixnormal);
	if (b314) {
		v[315] = 0e0;
		v[316] = 0e0;
		v[317] = 0e0;
		v[318] = 0e0;
		v[319] = 0e0;
		v[320] = 0e0;
		v[321] = 0e0;
		v[322] = 0e0;
		v[323] = 0.5e0*v[423];
		v[324] = 0.5e0*v[424];
		v[325] = 0.5e0*v[425];
	}
	else {
		v[327] = -0.5e0*v[422];
		v[326] = 0.5e0*v[421];
		v[320] = 0.5e0*i426*(v[296] * v[308] + v[299] * v[309] + v[301] * v[310]);
		v[319] = v[308] * v[326];
		v[316] = -0.5e0*i427*(v[289] * v[308] + v[292] * v[309] + v[294] * v[310]);
		v[315] = v[308] * v[327];
		v[323] = 0.5e0*v[428];
		v[324] = 0.5e0*v[429];
		v[325] = 0.5e0*v[430];
		v[321] = v[309] * v[326];
		v[322] = 1e0*v[310] * v[326];
		v[317] = v[309] * v[327];
		v[318] = 1e0*v[310] * v[327];
	};
	v[331] = -((v[291] * v[316]) / v[330]);
	v[329] = -((v[298] * v[320]) / v[328]);
	v[319] = v[319] + v[296] * v[329];
	v[321] = v[321] + v[299] * v[329];
	v[322] = v[322] + v[301] * v[329];
	v[315] = v[315] + v[289] * v[331];
	v[317] = v[317] + v[292] * v[331];
	v[318] = v[318] + v[294] * v[331];
	v[345] = (-(v[276] * v[319]) + v[275] * v[321]) / 2e0;
	v[333] = -(v[283] * v[319]) + v[280] * v[321];
	v[339] = (-(v[277] * v[321]) + v[276] * v[322]) / 2e0;
	v[342] = (v[277] * v[319] - v[275] * v[322]) / 2e0;
	v[336] = -(v[286] * v[321]) + v[283] * v[322];
	v[337] = v[286] * v[319] - v[280] * v[322];
	v[338] = -(v[236] * v[325]) + v[339];
	v[340] = -(v[237] * v[325]) - v[339];
	v[341] = -(v[236] * v[324]) + v[342];
	v[343] = -(v[237] * v[324]) - v[342];
	v[344] = -(v[236] * v[323]) + v[345];
	v[346] = -(v[237] * v[323]) - v[345];
	v[347] = v[213] * v[338] + v[222] * v[340] + v[216] * v[341] + v[225] * v[343] + v[219] * v[344] + v[228] * v[346];
	v[348] = v[212] * v[338] + v[221] * v[340] + v[215] * v[341] + v[224] * v[343] + v[218] * v[344] + v[227] * v[346];
	v[349] = (v[218] * v[236] + v[227] * v[237])*v[333] + (v[212] * v[236] + v[221] * v[237])*v[336] + (v[215] * v[236]
		+ v[224] * v[237])*v[337];
	v[438] = (*aB)*v[349];
	v[350] = (v[219] * v[236] + v[228] * v[237])*v[333] + (v[213] * v[236] + v[222] * v[237])*v[336] + (v[216] * v[236]
		+ v[225] * v[237])*v[337];
	v[439] = (*bB)*v[350];
	v[351] = v[267] * v[438] + v[242] * v[439];
	v[440] = ((*aB)*(v[267] * v[348] - v[242] * v[349]) + (*bB)*(v[242] * v[347] + v[267] * v[350]))*v[352] + Power
	(v[268], -2e0 - 1e0 / v[238])*v[351] * v[353] * v[355];
	v[431] = v[351] * v[352];
	v[406] = v[404] * v[431];
	v[402] = -(v[401] * v[431]);
	v[370] = (-(v[256] * v[315]) + v[255] * v[317]) / 2e0;
	v[358] = -(v[263] * v[315]) + v[260] * v[317];
	v[364] = (-(v[257] * v[317]) + v[256] * v[318]) / 2e0;
	v[367] = (v[257] * v[315] - v[255] * v[318]) / 2e0;
	v[361] = -(v[266] * v[317]) + v[263] * v[318];
	v[362] = v[266] * v[315] - v[260] * v[318];
	v[363] = v[163] * v[325] + v[364];
	v[365] = v[164] * v[325] - v[364];
	v[366] = v[163] * v[324] + v[367];
	v[368] = v[164] * v[324] - v[367];
	v[369] = v[163] * v[323] + v[370];
	v[371] = v[164] * v[323] - v[370];
	v[372] = v[140] * v[363] + v[149] * v[365] + v[143] * v[366] + v[152] * v[368] + v[146] * v[369] + v[155] * v[371];
	v[373] = v[139] * v[363] + v[148] * v[365] + v[142] * v[366] + v[151] * v[368] + v[145] * v[369] + v[154] * v[371];
	v[374] = (v[145] * v[163] + v[154] * v[164])*v[358] + (v[139] * v[163] + v[148] * v[164])*v[361] + (v[142] * v[163]
		+ v[151] * v[164])*v[362];
	v[435] = (*aA)*v[374];
	v[375] = (v[146] * v[163] + v[155] * v[164])*v[358] + (v[140] * v[163] + v[149] * v[164])*v[361] + (v[143] * v[163]
		+ v[152] * v[164])*v[362];
	v[436] = (*bA)*v[375];
	v[376] = v[247] * v[435] + v[169] * v[436];
	v[437] = ((*aA)*(v[247] * v[373] - v[169] * v[374]) + (*bA)*(v[169] * v[372] + v[247] * v[375]))*v[377] + Power
	(v[248], -2e0 - 1e0 / v[165])*v[376] * v[378] * v[380];
	v[433] = v[376] * v[377];
	v[393] = v[391] * v[433];
	v[389] = -(v[388] * v[433]);
	v[489] = (-(v[264] * v[323]) + v[265] * v[323] - v[261] * v[324] + v[262] * v[324] - v[258] * v[325] + v[259] * v[325]
		- v[358] * v[382] - v[361] * v[383] - v[362] * v[384] + v[358] * v[385] + v[361] * v[386] + v[362] * v[387]) / 2e0;
	v[490] = v[247] * (v[392] * v[393] + v[166] * ((*bA)*v[372] - v[435]) + v[252] * v[436] + v[388] * (v[247] * v[250] * Power
	(v[251], v[394])*v[389] - v[390] * v[437])) - v[169] * (v[389] * v[390] + v[252] * v[435] + v[166] * ((*aA)*v[373]
		+ v[436]) + v[391] * (v[169] * Power(v[249], v[394])*v[250] * v[393] - v[392] * v[437]));
	v[491] = (v[284] * v[323] - v[285] * v[323] + v[281] * v[324] - v[282] * v[324] + v[278] * v[325] - v[279] * v[325]
		- v[333] * v[395] - v[336] * v[396] - v[337] * v[397] + v[333] * v[398] + v[336] * v[399] + v[337] * v[400]) / 2e0;
	v[492] = v[267] * (v[405] * v[406] + v[239] * ((*bB)*v[347] - v[438]) + v[272] * v[439] + v[401] * (v[267] * v[270] * Power
	(v[271], v[407])*v[402] - v[403] * v[440])) - v[242] * (v[402] * v[403] + v[272] * v[438] + v[239] * ((*aB)*v[348]
		+ v[439]) + v[404] * (v[242] * Power(v[269], v[407])*v[270] * v[406] - v[405] * v[440]));
	
	for (i312 = 1; i312 <= 4; i312++) {
		Gra[i312 - 1] = v[488 + i312];
	};/* end for */


	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
		delete[] QBBi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
	delete[] QBBi;
}

//Calcula a Hessiana do gap
void FlexibleSECylinder_1_FlexibleSECylinder_1::HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variável temporária - AceGen
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
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
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	aA = &surf1->a;
	aB = &surf2->a;
	bA = &surf1->b;
	bB = &surf2->b;
	eA = &surf1->e;
	eB = &surf2->e;
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_A->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_A->getMatrix();
	dduiB = surf2->ddui_A->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_AAi->getMatrix();
	xBBi = surf2->x_BAi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	QBBi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
		QBBi[i] = new double[3];
	}
	//Salvando variáveis locais para montagem de superfícies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_AAi->MatrixToPtr(QABi, 3);
	surf2->Q_BAi->MatrixToPtr(QBBi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	double *c = mc.getMatrix();
	double Hes[4][4];
	//Pontos de singularidade - função phi
	double test1 = 2 * c[1] / PI;
	test1 = abs(test1 - int(test1));//testa o quão próximo de inteiro é esse número
	if (test1 < tol_small_1)
		c[1] += tol_small_1;
	double test3 = 2 * c[3] / PI;
	test3 = abs(test3 - int(test3));//testa o quão próximo de inteiro é esse número
	if (test3 < tol_small_1)
		c[3] += tol_small_1;
	
	
	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	//AceGen
	int i312, i411, i681, i684, i687, i690, i718, i732, b303, b314, b509, b641;
	v[247] = cos(c[1]);
	v[702] = (*aA)*v[247];
	v[169] = sin(c[1]);
	v[703] = (*bA)*v[169];
	v[610] = (v[169] * v[169]);
	v[388] = _copysign(1.e0, v[169]);
	v[251] = fabs(v[169]);
	v[391] = _copysign(1.e0, v[247]);
	v[249] = fabs(v[247]);
	v[267] = cos(c[3]);
	v[695] = (*aB)*v[267];
	v[242] = sin(c[3]);
	v[696] = (*bB)*v[242];
	v[627] = (v[242] * v[242]);
	v[401] = _copysign(1.e0, v[242]);
	v[271] = fabs(v[242]);
	v[404] = _copysign(1.e0, v[267]);
	v[269] = fabs(v[267]);
	v[109] = Power(dA[3], 2);
	v[107] = (dA[3] * dA[4]) / 2e0;
	v[102] = Power(dA[4], 2);
	v[114] = (dA[4] * dA[5]) / 2e0;
	v[112] = (dA[3] * dA[5]) / 2e0;
	v[103] = Power(dA[5], 2);
	v[643] = v[102] + v[103];
	v[128] = Power(dA[9], 2);
	v[126] = (dA[10] * dA[9]) / 2e0;
	v[121] = Power(dA[10], 2);
	v[133] = (dA[10] * dA[11]) / 2e0;
	v[131] = (dA[11] * dA[9]) / 2e0;
	v[122] = Power(dA[11], 2);
	v[644] = v[121] + v[122];
	v[182] = Power(dB[3], 2);
	v[180] = (dB[3] * dB[4]) / 2e0;
	v[175] = Power(dB[4], 2);
	v[187] = (dB[4] * dB[5]) / 2e0;
	v[185] = (dB[3] * dB[5]) / 2e0;
	v[176] = Power(dB[5], 2);
	v[650] = v[175] + v[176];
	v[201] = Power(dB[9], 2);
	v[199] = (dB[10] * dB[9]) / 2e0;
	v[194] = Power(dB[10], 2);
	v[206] = (dB[10] * dB[11]) / 2e0;
	v[204] = (dB[11] * dB[9]) / 2e0;
	v[195] = Power(dB[11], 2);
	v[651] = v[194] + v[195];
	v[101] = 4e0 / (4e0 + v[109] + v[643]);
	v[104] = 1e0 - (v[101] * v[643]) / 2e0;
	v[105] = v[101] * (-dA[5] + v[107]);
	v[106] = v[101] * (dA[4] + v[112]);
	v[108] = v[101] * (dA[5] + v[107]);
	v[110] = 1e0 - (v[101] * (v[103] + v[109])) / 2e0;
	v[111] = v[101] * (-dA[3] + v[114]);
	v[113] = v[101] * (-dA[4] + v[112]);
	v[115] = v[101] * (dA[3] + v[114]);
	v[116] = 1e0 - (v[101] * (v[102] + v[109])) / 2e0;
	v[120] = 4e0 / (4e0 + v[128] + v[644]);
	v[123] = 1e0 - (v[120] * v[644]) / 2e0;
	v[124] = v[120] * (-dA[11] + v[126]);
	v[125] = v[120] * (dA[10] + v[131]);
	v[127] = v[120] * (dA[11] + v[126]);
	v[129] = 1e0 - (v[120] * (v[122] + v[128])) / 2e0;
	v[130] = v[120] * (-dA[9] + v[133]);
	v[132] = v[120] * (-dA[10] + v[131]);
	v[134] = v[120] * (dA[9] + v[133]);
	v[135] = 1e0 - (v[120] * (v[121] + v[128])) / 2e0;
	v[139] = QAAi[0][0] * v[104] + QAAi[1][0] * v[105] + QAAi[2][0] * v[106];
	v[140] = QAAi[0][1] * v[104] + QAAi[1][1] * v[105] + QAAi[2][1] * v[106];
	v[142] = QAAi[0][0] * v[108] + QAAi[1][0] * v[110] + QAAi[2][0] * v[111];
	v[143] = QAAi[0][1] * v[108] + QAAi[1][1] * v[110] + QAAi[2][1] * v[111];
	v[145] = QAAi[0][0] * v[113] + QAAi[1][0] * v[115] + QAAi[2][0] * v[116];
	v[146] = QAAi[0][1] * v[113] + QAAi[1][1] * v[115] + QAAi[2][1] * v[116];
	v[148] = QBAi[0][0] * v[123] + QBAi[1][0] * v[124] + QBAi[2][0] * v[125];
	v[149] = QBAi[0][1] * v[123] + QBAi[1][1] * v[124] + QBAi[2][1] * v[125];
	v[151] = QBAi[0][0] * v[127] + QBAi[1][0] * v[129] + QBAi[2][0] * v[130];
	v[152] = QBAi[0][1] * v[127] + QBAi[1][1] * v[129] + QBAi[2][1] * v[130];
	v[154] = QBAi[0][0] * v[132] + QBAi[1][0] * v[134] + QBAi[2][0] * v[135];
	v[155] = QBAi[0][1] * v[132] + QBAi[1][1] * v[134] + QBAi[2][1] * v[135];
	v[163] = (1e0 - c[0]) / 2e0;
	v[164] = (1e0 + c[0]) / 2e0;
	v[437] = v[145] * v[163] + v[154] * v[164];
	v[436] = v[146] * v[163] + v[155] * v[164];
	v[434] = v[139] * v[163] + v[148] * v[164];
	v[433] = v[140] * v[163] + v[149] * v[164];
	v[430] = v[142] * v[163] + v[151] * v[164];
	v[429] = v[143] * v[163] + v[152] * v[164];
	v[165] = 2e0 / (*eA);
	v[613] = -(v[165] * v[169] * v[391]);
	v[607] = v[165] * v[247] * v[388];
	v[378] = -1e0 - 1e0 / v[165];
	v[596] = -1e0 + v[378];
	v[250] = -1e0 + v[165];
	v[605] = (v[247] * v[247])*v[250] * v[388];
	v[394] = -1e0 + v[250];
	v[611] = -1e0 + v[394];
	v[423] = Power(v[249], v[394]);
	v[679] = v[250] * v[423];
	v[609] = -(v[391] * v[679]);
	v[422] = Power(v[251], v[394]);
	v[719] = v[250] * v[422];
	v[392] = Power(v[249], v[250]);
	v[737] = v[247] * v[392] + v[609] * v[610];
	v[721] = v[391] * v[392];
	v[390] = Power(v[251], v[250]);
	v[738] = v[169] * v[390] - v[422] * v[605];
	v[722] = v[388] * v[390];
	v[380] = v[390] * v[607] + v[392] * v[613];
	v[597] = (v[380] * v[380]);
	v[248] = Power(v[249], v[165]) + Power(v[251], v[165]);
	v[756] = Power(v[248], -3e0 - 1e0 / v[165]);
	v[707] = Power(v[248], v[596])*v[378];
	v[421] = -(v[707] / v[165]);
	v[377] = Power(v[248], v[378]);
	v[678] = -(v[377] / v[165]);
	v[252] = v[380] * v[678];
	v[646] = (*aA)*v[252];
	v[645] = (*bA)*v[252];
	v[166] = 1e0 / Power(v[248], 1 / v[165]);
	v[649] = (*bA)*v[166];
	v[648] = (*aA)*v[166];
	v[254] = v[169] * v[645] + v[247] * v[649];
	v[253] = v[247] * v[646] - v[169] * v[648];
	v[387] = v[151] * v[253] + v[152] * v[254];
	v[386] = v[148] * v[253] + v[149] * v[254];
	v[385] = v[154] * v[253] + v[155] * v[254];
	v[384] = v[142] * v[253] + v[143] * v[254];
	v[717] = v[384] - v[387];
	v[383] = v[139] * v[253] + v[140] * v[254];
	v[716] = v[383] - v[386];
	v[382] = v[145] * v[253] + v[146] * v[254];
	v[715] = v[382] - v[385];
	v[257] = v[163] * v[382] + v[164] * v[385];
	v[256] = v[163] * v[384] + v[164] * v[387];
	v[255] = v[163] * v[383] + v[164] * v[386];
	v[168] = v[247] * v[648];
	v[740] = v[168] + 2e0*v[169] * v[646];
	v[170] = v[169] * v[649];
	v[739] = v[170] - 2e0*v[247] * v[645];
	v[265] = dA[8] + v[154] * v[168] + v[155] * v[170] + xBAi[2];
	v[264] = dA[2] + v[145] * v[168] + v[146] * v[170] + xAAi[2];
	v[714] = -v[264] + v[265];
	v[266] = v[714] / 2e0;
	v[262] = dA[7] + v[151] * v[168] + v[152] * v[170] + xBAi[1];
	v[261] = dA[1] + v[142] * v[168] + v[143] * v[170] + xAAi[1];
	v[713] = -v[261] + v[262];
	v[263] = v[713] / 2e0;
	v[259] = dA[6] + v[148] * v[168] + v[149] * v[170] + xBAi[0];
	v[258] = dA[0] + v[139] * v[168] + v[140] * v[170] + xAAi[0];
	v[712] = -v[258] + v[259];
	v[260] = v[712] / 2e0;
	v[294] = -(v[256] * v[260]) + v[255] * v[263];
	v[292] = v[257] * v[260] - v[255] * v[266];
	v[174] = 4e0 / (4e0 + v[182] + v[650]);
	v[177] = 1e0 - (v[174] * v[650]) / 2e0;
	v[178] = v[174] * (-dB[5] + v[180]);
	v[179] = v[174] * (dB[4] + v[185]);
	v[181] = v[174] * (dB[5] + v[180]);
	v[183] = 1e0 - (v[174] * (v[176] + v[182])) / 2e0;
	v[184] = v[174] * (-dB[3] + v[187]);
	v[186] = v[174] * (-dB[4] + v[185]);
	v[188] = v[174] * (dB[3] + v[187]);
	v[189] = 1e0 - (v[174] * (v[175] + v[182])) / 2e0;
	v[193] = 4e0 / (4e0 + v[201] + v[651]);
	v[196] = 1e0 - (v[193] * v[651]) / 2e0;
	v[197] = v[193] * (-dB[11] + v[199]);
	v[198] = v[193] * (dB[10] + v[204]);
	v[200] = v[193] * (dB[11] + v[199]);
	v[202] = 1e0 - (v[193] * (v[195] + v[201])) / 2e0;
	v[203] = v[193] * (-dB[9] + v[206]);
	v[205] = v[193] * (-dB[10] + v[204]);
	v[207] = v[193] * (dB[9] + v[206]);
	v[208] = 1e0 - (v[193] * (v[194] + v[201])) / 2e0;
	v[212] = QABi[0][0] * v[177] + QABi[1][0] * v[178] + QABi[2][0] * v[179];
	v[213] = QABi[0][1] * v[177] + QABi[1][1] * v[178] + QABi[2][1] * v[179];
	v[215] = QABi[0][0] * v[181] + QABi[1][0] * v[183] + QABi[2][0] * v[184];
	v[216] = QABi[0][1] * v[181] + QABi[1][1] * v[183] + QABi[2][1] * v[184];
	v[218] = QABi[0][0] * v[186] + QABi[1][0] * v[188] + QABi[2][0] * v[189];
	v[219] = QABi[0][1] * v[186] + QABi[1][1] * v[188] + QABi[2][1] * v[189];
	v[221] = QBBi[0][0] * v[196] + QBBi[1][0] * v[197] + QBBi[2][0] * v[198];
	v[222] = QBBi[0][1] * v[196] + QBBi[1][1] * v[197] + QBBi[2][1] * v[198];
	v[224] = QBBi[0][0] * v[200] + QBBi[1][0] * v[202] + QBBi[2][0] * v[203];
	v[225] = QBBi[0][1] * v[200] + QBBi[1][1] * v[202] + QBBi[2][1] * v[203];
	v[227] = QBBi[0][0] * v[205] + QBBi[1][0] * v[207] + QBBi[2][0] * v[208];
	v[228] = QBBi[0][1] * v[205] + QBBi[1][1] * v[207] + QBBi[2][1] * v[208];
	v[236] = (1e0 - c[2]) / 2e0;
	v[237] = (1e0 + c[2]) / 2e0;
	v[471] = v[218] * v[236] + v[227] * v[237];
	v[470] = v[219] * v[236] + v[228] * v[237];
	v[468] = v[212] * v[236] + v[221] * v[237];
	v[467] = v[213] * v[236] + v[222] * v[237];
	v[464] = v[215] * v[236] + v[224] * v[237];
	v[463] = v[216] * v[236] + v[225] * v[237];
	v[238] = 2e0 / (*eB);
	v[630] = -(v[238] * v[242] * v[404]);
	v[624] = v[238] * v[267] * v[401];
	v[353] = -1e0 - 1e0 / v[238];
	v[553] = -1e0 + v[353];
	v[270] = -1e0 + v[238];
	v[622] = (v[267] * v[267])*v[270] * v[401];
	v[407] = -1e0 + v[270];
	v[628] = -1e0 + v[407];
	v[457] = Power(v[269], v[407]);
	v[671] = v[270] * v[457];
	v[626] = -(v[404] * v[671]);
	v[456] = Power(v[271], v[407]);
	v[733] = v[270] * v[456];
	v[405] = Power(v[269], v[270]);
	v[741] = v[267] * v[405] + v[626] * v[627];
	v[735] = v[404] * v[405];
	v[403] = Power(v[271], v[270]);
	v[742] = v[242] * v[403] - v[456] * v[622];
	v[736] = v[401] * v[403];
	v[355] = v[403] * v[624] + v[405] * v[630];
	v[554] = (v[355] * v[355]);
	v[268] = Power(v[269], v[238]) + Power(v[271], v[238]);
	v[754] = Power(v[268], -3e0 - 1e0 / v[238]);
	v[700] = Power(v[268], v[553])*v[353];
	v[455] = -(v[700] / v[238]);
	v[352] = Power(v[268], v[353]);
	v[670] = -(v[352] / v[238]);
	v[272] = v[355] * v[670];
	v[653] = (*aB)*v[272];
	v[652] = (*bB)*v[272];
	v[239] = 1e0 / Power(v[268], 1 / v[238]);
	v[656] = (*bB)*v[239];
	v[655] = (*aB)*v[239];
	v[274] = v[242] * v[652] + v[267] * v[656];
	v[273] = v[267] * v[653] - v[242] * v[655];
	v[400] = v[224] * v[273] + v[225] * v[274];
	v[399] = v[221] * v[273] + v[222] * v[274];
	v[398] = v[227] * v[273] + v[228] * v[274];
	v[397] = v[215] * v[273] + v[216] * v[274];
	v[731] = -v[397] + v[400];
	v[396] = v[212] * v[273] + v[213] * v[274];
	v[730] = v[396] - v[399];
	v[395] = v[218] * v[273] + v[219] * v[274];
	v[729] = v[395] - v[398];
	v[277] = v[236] * v[395] + v[237] * v[398];
	v[276] = v[236] * v[397] + v[237] * v[400];
	v[275] = v[236] * v[396] + v[237] * v[399];
	v[241] = v[267] * v[655];
	v[744] = v[241] + 2e0*v[242] * v[653];
	v[243] = v[242] * v[656];
	v[743] = v[243] - 2e0*v[267] * v[652];
	v[285] = dB[8] + v[227] * v[241] + v[228] * v[243] + xBBi[2];
	v[284] = dB[2] + v[218] * v[241] + v[219] * v[243] + xABi[2];
	v[728] = -v[284] + v[285];
	v[286] = v[728] / 2e0;
	v[282] = dB[7] + v[224] * v[241] + v[225] * v[243] + xBBi[1];
	v[281] = dB[1] + v[215] * v[241] + v[216] * v[243] + xABi[1];
	v[727] = -v[281] + v[282];
	v[283] = v[727] / 2e0;
	v[279] = dB[6] + v[221] * v[241] + v[222] * v[243] + xBBi[0];
	v[278] = dB[0] + v[212] * v[241] + v[213] * v[243] + xABi[0];
	v[726] = -v[278] + v[279];
	v[280] = v[726] / 2e0;
	v[301] = -(v[276] * v[280]) + v[275] * v[283];
	v[299] = v[277] * v[280] - v[275] * v[286];
	v[308] = v[163] * v[258] + v[164] * v[259] - v[236] * v[278] - v[237] * v[279];
	v[309] = v[163] * v[261] + v[164] * v[262] - v[236] * v[281] - v[237] * v[282];
	v[310] = v[163] * v[264] + v[164] * v[265] - v[236] * v[284] - v[237] * v[285];
	v[287] = ((*normalintA) ? -1 : 1);
	v[288] = ((*normalintB) ? -1 : 1);
	v[289] = -(v[257] * v[263]) + v[256] * v[266];
	v[330] = (v[289] * v[289]) + (v[292] * v[292]) + (v[294] * v[294]);
	v[752] = 1e0 / Power(v[330], 2);
	v[291] = 1e0 / sqrt(v[330]);
	v[296] = -(v[277] * v[283]) + v[276] * v[286];
	v[328] = (v[296] * v[296]) + (v[299] * v[299]) + (v[301] * v[301]);
	v[751] = 1e0 / Power(v[328], 2);
	v[298] = 1e0 / sqrt(v[328]);
	if ((*fixnormal)) {
		v[661] = -normalA[0] + normalB[0];
		v[660] = -normalA[1] + normalB[1];
		v[659] = -normalA[2] + normalB[2];
		v[307] = 0.5e0*(v[310] * v[659] + v[309] * v[660] + v[308] * v[661]);
	}
	else {
		v[933] = v[266];
		v[934] = 0e0;
		v[935] = -v[286];
		v[936] = 0e0;
		v[937] = v[263];
		v[938] = 0e0;
		v[939] = -v[283];
		v[940] = 0e0;
		v[941] = v[260];
		v[942] = 0e0;
		v[943] = -v[280];
		v[944] = 0e0;
		v[658] = v[287] * v[291];
		v[657] = v[288] * v[298];
		v[662] = v[301] * v[657] - v[294] * v[658];
		v[663] = v[299] * v[657] - v[292] * v[658];
		v[664] = v[296] * v[657] - v[289] * v[658];
		v[307] = 0.5e0*(v[310] * v[662] + v[309] * v[663] + v[308] * v[664]);
	};
	b314 = (*fixnormal);
	if (b314) {
		v[315] = 0e0;
		v[316] = 0e0;
		v[317] = 0e0;
		v[318] = 0e0;
		v[319] = 0e0;
		v[320] = 0e0;
		v[321] = 0e0;
		v[322] = 0e0;
		v[323] = 0.5e0*v[659];
		v[324] = 0.5e0*v[660];
		v[325] = 0.5e0*v[661];
	}
	else {
		v[327] = -0.5e0*v[658];
		v[326] = 0.5e0*v[657];
		v[320] = 0.5e0*v[288] * (v[296] * v[308] + v[299] * v[309] + v[301] * v[310]);
		v[319] = v[308] * v[326];
		v[316] = -0.5e0*v[287] * (v[289] * v[308] + v[292] * v[309] + v[294] * v[310]);
		v[315] = v[308] * v[327];
		v[323] = 0.5e0*v[662];
		v[324] = 0.5e0*v[663];
		v[325] = 0.5e0*v[664];
		v[321] = v[309] * v[326];
		v[322] = 1e0*v[310] * v[326];
		v[317] = v[309] * v[327];
		v[318] = 1e0*v[310] * v[327];
	};
	v[694] = -(v[316] / v[330]);
	v[331] = v[291] * v[694];
	v[693] = -(v[320] / v[328]);
	v[329] = v[298] * v[693];
	v[319] = v[319] + v[296] * v[329];
	v[321] = v[321] + v[299] * v[329];
	v[322] = v[322] + v[301] * v[329];
	v[315] = v[315] + v[289] * v[331];
	v[317] = v[317] + v[292] * v[331];
	v[318] = v[318] + v[294] * v[331];
	v[345] = (-(v[276] * v[319]) + v[275] * v[321]) / 2e0;
	v[333] = -(v[283] * v[319]) + v[280] * v[321];
	v[339] = (-(v[277] * v[321]) + v[276] * v[322]) / 2e0;
	v[342] = (v[277] * v[319] - v[275] * v[322]) / 2e0;
	v[336] = -(v[286] * v[321]) + v[283] * v[322];
	v[337] = v[286] * v[319] - v[280] * v[322];
	v[762] = (-v[218] + v[227])*v[333] + (-v[212] + v[221])*v[336] + (-v[215] + v[224])*v[337];
	v[761] = (-v[219] + v[228])*v[333] + (-v[213] + v[222])*v[336] + (-v[216] + v[225])*v[337];
	v[338] = -(v[236] * v[325]) + v[339];
	v[340] = -(v[237] * v[325]) - v[339];
	v[341] = -(v[236] * v[324]) + v[342];
	v[343] = -(v[237] * v[324]) - v[342];
	v[344] = -(v[236] * v[323]) + v[345];
	v[346] = -(v[237] * v[323]) - v[345];
	v[667] = (*bB)*(v[213] * v[338] + v[222] * v[340] + v[216] * v[341] + v[225] * v[343] + v[219] * v[344] + v[228] * v[346]);
	v[348] = v[212] * v[338] + v[221] * v[340] + v[215] * v[341] + v[224] * v[343] + v[218] * v[344] + v[227] * v[346];
	v[668] = (*aB)*v[348];
	v[632] = v[348] * v[655];
	v[349] = v[337] * v[464] + v[336] * v[468] + v[333] * v[471];
	v[665] = -((*aB)*v[349]);
	v[633] = v[349] * v[653];
	v[631] = v[665] + v[667];
	v[350] = v[337] * v[463] + v[336] * v[467] + v[333] * v[470];
	v[666] = (*bB)*v[350];
	v[634] = v[350] * v[656];
	v[620] = v[666] + v[668];
	v[753] = -(v[242] * v[620]) + v[267] * v[631];
	v[697] = v[242] * v[665] + v[267] * v[666];
	v[558] = v[242] * v[667] + v[267] * v[668] + v[697];
	v[351] = -(v[267] * v[665]) + v[242] * v[666];
	v[669] = v[351] * v[352];
	v[551] = v[351] * v[355] * v[455];
	v[406] = v[404] * v[669];
	v[734] = v[242] * v[406];
	v[402] = -(v[401] * v[669]);
	v[635] = v[402] * v[403];
	v[623] = v[267] * v[402] * v[733];
	v[356] = v[551] + v[558] * v[670];
	v[672] = v[238] * v[356];
	v[636] = v[404] * (v[405] * v[672] + v[671] * v[734]);
	v[764] = v[272] * v[631] - v[632] - v[633] - v[634] - v[635] - v[636] + v[401] * (Power(v[271], v[628]
	)*v[402] * v[407] * v[622] + (-(v[242] * v[402]) + v[356] * v[624])*v[733]) + v[626] * v[734];
	v[621] = v[405] * v[406] + v[239] * v[631] + v[272] * v[666] + v[401] * (v[623] + v[403] * v[672]);
	v[763] = v[272] * v[620] + v[621] + v[401] * v[623] + v[270] * v[404] * (-(Power(v[269], v[628]
	)*v[404] * v[406] * v[407] * v[627]) + v[457] * (v[267] * v[406] + v[356] * v[630]));
	v[370] = (-(v[256] * v[315]) + v[255] * v[317]) / 2e0;
	v[358] = -(v[263] * v[315]) + v[260] * v[317];
	v[364] = (-(v[257] * v[317]) + v[256] * v[318]) / 2e0;
	v[367] = (v[257] * v[315] - v[255] * v[318]) / 2e0;
	v[361] = -(v[266] * v[317]) + v[263] * v[318];
	v[362] = v[266] * v[315] - v[260] * v[318];
	v[758] = (-v[145] + v[154])*v[358] + (-v[139] + v[148])*v[361] + (-v[142] + v[151])*v[362];
	v[757] = (-v[146] + v[155])*v[358] + (-v[140] + v[149])*v[361] + (-v[143] + v[152])*v[362];
	v[363] = v[163] * v[325] + v[364];
	v[365] = v[164] * v[325] - v[364];
	v[366] = v[163] * v[324] + v[367];
	v[368] = v[164] * v[324] - v[367];
	v[369] = v[163] * v[323] + v[370];
	v[371] = v[164] * v[323] - v[370];
	v[675] = (*bA)*(v[140] * v[363] + v[149] * v[365] + v[143] * v[366] + v[152] * v[368] + v[146] * v[369] + v[155] * v[371]);
	v[373] = v[139] * v[363] + v[148] * v[365] + v[142] * v[366] + v[151] * v[368] + v[145] * v[369] + v[154] * v[371];
	v[676] = (*aA)*v[373];
	v[615] = v[373] * v[648];
	v[374] = v[362] * v[430] + v[361] * v[434] + v[358] * v[437];
	v[673] = -((*aA)*v[374]);
	v[616] = v[374] * v[646];
	v[614] = v[673] + v[675];
	v[375] = v[362] * v[429] + v[361] * v[433] + v[358] * v[436];
	v[674] = (*bA)*v[375];
	v[617] = v[375] * v[649];
	v[603] = v[674] + v[676];
	v[755] = -(v[169] * v[603]) + v[247] * v[614];
	v[704] = v[169] * v[673] + v[247] * v[674];
	v[601] = v[169] * v[675] + v[247] * v[676] + v[704];
	v[376] = -(v[247] * v[673]) + v[169] * v[674];
	v[677] = v[376] * v[377];
	v[594] = v[376] * v[380] * v[421];
	v[393] = v[391] * v[677];
	v[720] = v[169] * v[393];
	v[389] = -(v[388] * v[677]);
	v[618] = v[389] * v[390];
	v[606] = v[247] * v[389] * v[719];
	v[381] = v[594] + v[601] * v[678];
	v[680] = v[165] * v[381];
	v[619] = v[391] * (v[392] * v[680] + v[679] * v[720]);
	v[760] = v[252] * v[614] - v[615] - v[616] - v[617] - v[618] - v[619] + v[388] * (Power(v[251], v[611]
	)*v[389] * v[394] * v[605] + (-(v[169] * v[389]) + v[381] * v[607])*v[719]) + v[609] * v[720];
	v[604] = v[392] * v[393] + v[166] * v[614] + v[252] * v[674] + v[388] * (v[606] + v[390] * v[680]);
	v[759] = v[252] * v[603] + v[604] + v[388] * v[606] + v[250] * v[391] * (-(Power(v[249], v[611]
	)*v[391] * v[393] * v[394] * v[610]) + v[423] * (v[247] * v[393] + v[381] * v[613]));
	v[841] = (-(v[264] * v[323]) + v[265] * v[323] - v[261] * v[324] + v[262] * v[324] - v[258] * v[325] + v[259] * v[325]
		- v[358] * v[382] - v[361] * v[383] - v[362] * v[384] + v[358] * v[385] + v[361] * v[386] + v[362] * v[387]) / 2e0;
	v[842] = v[247] * v[604] - v[169] * (v[615] + v[616] + v[617] + v[618] + v[619]);
	v[843] = (v[284] * v[323] - v[285] * v[323] + v[281] * v[324] - v[282] * v[324] + v[278] * v[325] - v[279] * v[325]
		- v[333] * v[395] - v[336] * v[396] - v[337] * v[397] + v[333] * v[398] + v[336] * v[399] + v[337] * v[400]) / 2e0;
	v[844] = v[267] * v[621] - v[242] * (v[632] + v[633] + v[634] + v[635] + v[636]);
	
	for (i312 = 1; i312 <= 4; i312++) {
		i732 = i312 == 4;
		i718 = i312 == 2;
		i690 = (i732 ? 1 : 0);
		v[698] = i690 * v[554];
		i687 = (i718 ? 1 : 0);
		v[705] = i687 * v[597];
		i684 = (i312 == 1 ? 1 : 0);
		v[689] = -i684 / 2e0;
		v[688] = -i684 / 2e0;
		v[686] = -i684 / 2e0;
		v[685] = -i684 / 2e0;
		i681 = (i312 == 3 ? 1 : 0);
		v[692] = -i681 / 2e0;
		v[691] = i681 / 2e0;
		v[683] = i681 / 2e0;
		v[682] = -i681 / 2e0;
		v[544] = -(i681*v[333]) / 2e0;
		v[541] = v[337] * v[682];
		v[538] = v[336] * v[682];
		v[532] = v[323] * v[683];
		v[528] = v[324] * v[683];
		v[524] = v[325] * v[691];
		v[587] = -(i684*v[358]) / 2e0;
		v[584] = v[362] * v[685];
		v[581] = v[361] * v[685];
		v[575] = v[323] * v[686];
		v[571] = v[324] * v[686];
		v[567] = v[325] * v[688];
		v[708] = i687 * (v[391] * v[737] + v[388] * v[738]);
		v[417] = i687 * v[253];
		v[420] = i687 * v[254];
		v[424] = v[421] * v[705] + v[377] * v[708];
		v[426] = v[424] * v[703] - i687 * v[739];
		v[428] = v[424] * v[702] - i687 * v[740];
		v[431] = v[426] * v[429] + v[428] * v[430] + v[688] * v[717];
		v[435] = v[426] * v[433] + v[428] * v[434] + v[689] * v[716];
		v[438] = v[426] * v[436] + v[428] * v[437] + v[689] * v[715];
		v[439] = v[154] * v[417] + v[155] * v[420];
		v[440] = v[145] * v[417] + v[146] * v[420];
		v[709] = v[439] - v[440];
		v[441] = v[151] * v[417] + v[152] * v[420];
		v[442] = v[142] * v[417] + v[143] * v[420];
		v[710] = v[441] - v[442];
		v[443] = v[148] * v[417] + v[149] * v[420];
		v[444] = v[139] * v[417] + v[140] * v[420];
		v[711] = v[443] - v[444];
		v[445] = v[710] / 2e0;
		v[446] = v[711] / 2e0;
		v[486] = -(v[260] * v[431]) + v[263] * v[435] + v[255] * v[445] - v[256] * v[446];
		v[447] = v[709] / 2e0;
		v[492] = v[266] * v[431] - v[263] * v[438] - v[257] * v[445] + v[256] * v[447];
		v[489] = -(v[266] * v[435]) + v[260] * v[438] + v[257] * v[446] - v[255] * v[447];
		v[701] = i690 * (v[404] * v[741] + v[401] * v[742]);
		v[451] = i690 * v[273];
		v[454] = i690 * v[274];
		v[458] = v[455] * v[698] + v[352] * v[701];
		v[460] = v[458] * v[696] - i690 * v[743];
		v[462] = v[458] * v[695] - i690 * v[744];
		v[465] = v[460] * v[463] + v[462] * v[464] + v[691] * v[731];
		v[469] = v[460] * v[467] + v[462] * v[468] + v[692] * v[730];
		v[472] = v[460] * v[470] + v[462] * v[471] + v[692] * v[729];
		v[473] = v[227] * v[451] + v[228] * v[454];
		v[474] = v[218] * v[451] + v[219] * v[454];
		v[723] = v[473] - v[474];
		v[475] = v[224] * v[451] + v[225] * v[454];
		v[476] = v[215] * v[451] + v[216] * v[454];
		v[724] = v[475] - v[476];
		v[477] = v[221] * v[451] + v[222] * v[454];
		v[478] = v[212] * v[451] + v[213] * v[454];
		v[725] = v[477] - v[478];
		v[482] = v[724] / 2e0;
		v[483] = v[725] / 2e0;
		v[496] = -(v[280] * v[465]) + v[283] * v[469] + v[275] * v[482] - v[276] * v[483];
		v[484] = v[723] / 2e0;
		v[502] = v[286] * v[465] - v[283] * v[472] - v[277] * v[482] + v[276] * v[484];
		v[499] = -(v[286] * v[469]) + v[280] * v[472] + v[277] * v[483] - v[275] * v[484];
		v[485] = v[331] * v[486];
		v[488] = v[331] * v[489];
		v[491] = v[294] * v[486] + v[292] * v[489] + v[289] * v[492];
		v[493] = v[331] * v[492];
		v[495] = v[329] * v[496];
		v[498] = v[329] * v[499];
		v[501] = v[301] * v[496] + v[299] * v[499] + v[296] * v[502];
		v[503] = v[329] * v[502];
		v[506] = v[501] * v[693];
		v[508] = v[491] * v[694];
		b509 = (*fixnormal);
		if (b509) {
			v[510] = 0e0;
			v[511] = 0e0;
			v[512] = 0e0;
		}
		else {
			v[481] = v[164] * v[443] + v[163] * v[444] - v[237] * v[477] - v[236] * v[478] + v[940 + i312];
			v[480] = v[164] * v[441] + v[163] * v[442] - v[237] * v[475] - v[236] * v[476] + v[936 + i312];
			v[479] = v[164] * v[439] + v[163] * v[440] - v[237] * v[473] - v[236] * v[474] + v[932 + i312];
			v[514] = (-0.5e0*v[501] * v[657]) / v[328];
			v[513] = (0.5e0*v[491] * v[658]) / v[330];
			v[493] = v[327] * v[481] + v[493] + v[308] * v[513];
			v[488] = v[327] * v[480] + v[488] + v[309] * v[513];
			v[485] = v[327] * v[479] + v[485] + 1e0*v[310] * v[513];
			v[503] = v[326] * v[481] + v[503] + v[308] * v[514];
			v[512] = 1e0*v[327] * v[486] + 1e0*v[326] * v[496] + 1e0*v[294] * v[513] + v[301] * v[514];
			v[511] = v[327] * v[489] + v[326] * v[499] + 1e0*v[292] * v[513] + 1e0*v[299] * v[514];
			v[510] = v[327] * v[492] + v[326] * v[502] + 1e0*v[289] * v[513] + 1e0*v[296] * v[514];
			v[498] = v[326] * v[480] + v[498] + 1e0*v[309] * v[514];
			v[495] = v[326] * v[479] + v[495] + 1e0*v[310] * v[514];
			v[506] = v[288] * (0.5e0*v[301] * v[479] + 0.5e0*v[299] * v[480] + 0.5e0*v[296] * v[481] + 0.5e0*v[310] * v[496]
				+ 0.5e0*v[309] * v[499] + 0.5e0*v[308] * v[502]) + v[506];
			v[508] = v[287] * (-0.5e0*v[294] * v[479] - 0.5e0*v[292] * v[480] - 0.5e0*v[289] * v[481] - 0.5e0*v[310] * v[486]
				- 0.5e0*v[309] * v[489] - 0.5e0*v[308] * v[492]) + v[508];
		};
		v[515] = -(v[298] * (-2e0*v[320] * v[501] + v[328] * v[506])*v[751]) / 2e0;
		v[503] = v[503] + 2e0*v[296] * v[515];
		v[498] = v[498] + 2e0*v[299] * v[515];
		v[495] = v[495] + 2e0*v[301] * v[515];
		v[516] = -(v[291] * (-2e0*v[316] * v[491] + v[330] * v[508])*v[752]) / 2e0;
		v[493] = v[493] + 2e0*v[289] * v[516];
		v[488] = v[488] + 2e0*v[292] * v[516];
		v[485] = v[485] + 2e0*v[294] * v[516];
		v[533] = (-(v[319] * v[465]) + v[321] * v[469] + v[275] * v[498] - v[276] * v[503]) / 2e0;
		v[518] = -(v[319] * v[482]) + v[321] * v[483] + v[280] * v[498] - v[283] * v[503];
		v[525] = (v[322] * v[465] - v[321] * v[472] + v[276] * v[495] - v[277] * v[498]) / 2e0;
		v[529] = (-(v[322] * v[469]) + v[319] * v[472] - v[275] * v[495] + v[277] * v[503]) / 2e0;
		v[521] = v[322] * v[482] - v[321] * v[484] + v[283] * v[495] - v[286] * v[498];
		v[522] = -(v[322] * v[483]) + v[319] * v[484] - v[280] * v[495] + v[286] * v[503];
		v[523] = -(v[236] * v[510]) + v[524] + v[525];
		v[526] = -(v[237] * v[510]) - v[524] - v[525];
		v[527] = -(v[236] * v[511]) + v[528] + v[529];
		v[530] = -(v[237] * v[511]) - v[528] - v[529];
		v[531] = -(v[236] * v[512]) + v[532] + v[533];
		v[534] = -(v[237] * v[512]) - v[532] - v[533];
		v[535] = v[213] * v[523] + v[222] * v[526] + v[216] * v[527] + v[225] * v[530] + v[219] * v[531] + v[228] * v[534];
		v[536] = v[212] * v[523] + v[221] * v[526] + v[215] * v[527] + v[224] * v[530] + v[218] * v[531] + v[227] * v[534];
		v[537] = v[236] * v[521] + v[538];
		v[539] = v[237] * v[521] - v[538];
		v[540] = v[236] * v[522] + v[541];
		v[542] = v[237] * v[522] - v[541];
		v[543] = v[236] * v[518] + v[544];
		v[545] = v[237] * v[518] - v[544];
		v[546] = v[212] * v[537] + v[221] * v[539] + v[215] * v[540] + v[224] * v[542] + v[218] * v[543] + v[227] * v[545];
		v[547] = v[213] * v[537] + v[222] * v[539] + v[216] * v[540] + v[225] * v[542] + v[219] * v[543] + v[228] * v[545];
		v[550] = v[546] * v[695] + v[547] * v[696] + i690 * v[697];
		v[552] = i690 * v[551] + v[550] * v[670];
		v[559] = -((v[700] * (v[355] * v[550] + i690 * v[355] * v[558] - v[238] * v[351] * v[701]) + v[352] * ((*aB)*
			(v[267] * v[536] - v[242] * v[546]) + (*bB)*(v[242] * v[535] + v[267] * v[547]) + i690 * v[753])
			+ v[351] * v[353] * v[553] * v[698] * v[754]) / v[238]);
		v[576] = (-(v[315] * v[431]) + v[317] * v[435] + v[255] * v[488] - v[256] * v[493]) / 2e0;
		v[561] = -(v[315] * v[445]) + v[317] * v[446] + v[260] * v[488] - v[263] * v[493];
		v[568] = (v[318] * v[431] - v[317] * v[438] + v[256] * v[485] - v[257] * v[488]) / 2e0;
		v[572] = (-(v[318] * v[435]) + v[315] * v[438] - v[255] * v[485] + v[257] * v[493]) / 2e0;
		v[564] = v[318] * v[445] - v[317] * v[447] + v[263] * v[485] - v[266] * v[488];
		v[565] = -(v[318] * v[446]) + v[315] * v[447] - v[260] * v[485] + v[266] * v[493];
		v[566] = v[163] * v[510] + v[567] + v[568];
		v[569] = v[164] * v[510] - v[567] - v[568];
		v[570] = v[163] * v[511] + v[571] + v[572];
		v[573] = v[164] * v[511] - v[571] - v[572];
		v[574] = v[163] * v[512] + v[575] + v[576];
		v[577] = v[164] * v[512] - v[575] - v[576];
		v[578] = v[140] * v[566] + v[149] * v[569] + v[143] * v[570] + v[152] * v[573] + v[146] * v[574] + v[155] * v[577];
		v[579] = v[139] * v[566] + v[148] * v[569] + v[142] * v[570] + v[151] * v[573] + v[145] * v[574] + v[154] * v[577];
		v[580] = v[163] * v[564] + v[581];
		v[582] = v[164] * v[564] - v[581];
		v[583] = v[163] * v[565] + v[584];
		v[585] = v[164] * v[565] - v[584];
		v[586] = v[163] * v[561] + v[587];
		v[588] = v[164] * v[561] - v[587];
		v[589] = v[139] * v[580] + v[148] * v[582] + v[142] * v[583] + v[151] * v[585] + v[145] * v[586] + v[154] * v[588];
		v[590] = v[140] * v[580] + v[149] * v[582] + v[143] * v[583] + v[152] * v[585] + v[146] * v[586] + v[155] * v[588];
		v[593] = v[589] * v[702] + v[590] * v[703] + i687 * v[704];
		v[595] = i687 * v[594] + v[593] * v[678];
		v[602] = -((v[707] * (v[380] * v[593] + i687 * v[380] * v[601] - v[165] * v[376] * v[708]) + v[377] * ((*aA)*
			(v[247] * v[579] - v[169] * v[589]) + (*bA)*(v[169] * v[578] + v[247] * v[590]) + i687 * v[755])
			+ v[376] * v[378] * v[596] * v[705] * v[756]) / v[165]);
		v[1121] = (v[323] * v[709] + v[324] * v[710] + v[325] * v[711] + v[510] * v[712] + v[511] * v[713] + v[512] * v[714]
			- v[561] * v[715] - v[564] * v[716] - v[565] * v[717] + v[426] * v[757] + v[428] * v[758]) / 2e0;
		v[1122] = v[247] * ((i718 ? v[760] : 0e0) + (*bA)*(v[375] * v[424] + v[166] * v[578] + v[252] * v[590]) - v[589] * v[648]
			+ v[595] * (v[388] * v[607] * v[719] - v[165] * v[721]) + v[165] * v[602] * v[722]) - v[169] * ((i718 ? v[759] : 0e0) + (*aA
				)*(v[374] * v[424] + v[166] * v[579] + v[252] * v[589]) + v[590] * v[649] + v[165] * v[602] * v[721] + v[595] * (-
				(v[609] * v[613]) + v[165] * v[722]));
		v[1123] = (-(v[323] * v[723]) - v[324] * v[724] - v[325] * v[725] - v[510] * v[726] - v[511] * v[727] - v[512] * v[728]
			- v[518] * v[729] - v[521] * v[730] + v[522] * v[731] + v[460] * v[761] + v[462] * v[762]) / 2e0;
		v[1124] = v[267] * ((i732 ? v[764] : 0e0) + (*bB)*(v[350] * v[458] + v[239] * v[535] + v[272] * v[547]) - v[546] * v[655]
			+ v[552] * (v[401] * v[624] * v[733] - v[238] * v[735]) + v[238] * v[559] * v[736]) - v[242] * ((i732 ? v[763] : 0e0) + (*aB
				)*(v[349] * v[458] + v[239] * v[536] + v[272] * v[546]) + v[547] * v[656] + v[238] * v[559] * v[735] + v[552] * (-
				(v[626] * v[630]) + v[238] * v[736]));
		
		for (i411 = i312; i411 <= 4; i411++) {
			v[638] = v[1120 + i411];
			Hes[i312 - 1][i411 - 1] = v[638];
			if (i312 != i411) {
				Hes[i411 - 1][i312 - 1] = v[638];
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
		delete[] QBBi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
	delete[] QBBi;
}

//Chute inicial para coordenadas convectivas do par de superfícies
void FlexibleSECylinder_1_FlexibleSECylinder_1::InitialGuess(SSContactData* c_data)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double tol_ortho = 1e-12;
	
	Matrix uAA(3);
	Matrix uAB(3);
	Matrix uBA(3);
	Matrix uBB(3);
	for (int i = 0; i < 3; i++)
	{
		uAA(i, 0) = (*surf1->d_A)(i, 0);
		uBA(i, 0) = (*surf1->d_A)(i + 6, 0);
		uAB(i, 0) = (*surf2->d_A)(i, 0);
		uBB(i, 0) = (*surf2->d_A)(i + 6, 0);
	}
	
	Matrix bA = 0.5*(*surf1->x_AAi + uAA + *surf1->x_BAi + uBA);
	Matrix tA = 0.5*(*surf1->x_BAi + uBA - *surf1->x_AAi - uAA);
	Matrix bB = 0.5*(*surf2->x_AAi + uAB + *surf2->x_BAi + uBB);
	Matrix tB = 0.5*(*surf2->x_BAi + uBB - *surf2->x_AAi - uAB);

	double csi_A, csi_B;
	//Se não houver paralelismo
	if (abs(dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)) > tol_ortho)
	{
		//Estimativa dos csi's com base na configuração atual - considerando-se que são elementos retilíneos (baseado em Wriggers e Zavarise, 1997)
		csi_A = dot(bA - bB, (1.0 / (dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)))*(tB*dot(tA, tB) - tA*dot(tB, tB)));
		csi_B = -1.0*dot(bA - bB, (1.0 / (dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)))*(tA*dot(tA, tB) - tB*dot(tA, tA)));
	}
	else//Caso haja paralelismo
	{
		Matrix temp_xAA = *surf1->x_AAi + uAA;
		Matrix temp_xBA = *surf1->x_BAi + uBA;
		Matrix temp_xAB = *surf2->x_AAi + uAB;
		Matrix temp_xBB = *surf2->x_BAi + uBB;
		bool inverted = false;
		//Se a orientação tangente relativa das barras é oposta, inverte os nós da barra B
		if (dot(tA, tB) < 0)
		{
			temp_xAB = *surf2->x_BAi + uBB;
			temp_xBB = *surf2->x_AAi + uAB;
			inverted = true;
		}
		//Versor tangente
		Matrix t = (1.0 / (norm(temp_xBB - temp_xAB)))*(temp_xBB - temp_xAB);
		//Comprimentos das barras
		double LA = norm(temp_xAA - temp_xBA);
		double LB = norm(temp_xAB - temp_xBB);
		double csiAi = dot(temp_xAB - bA, t)*2.0 / LA;
		double csiAf = dot(temp_xBB - bA, t)*2.0 / LA;
		double csiBi = dot(temp_xAA - bB, t)*2.0 / LB;
		double csiBf = dot(temp_xBA - bB, t)*2.0 / LB;
		//Verificação de fora de range:
		if (csiAf > +1)
			csiAf = +1;
		if (csiAf < -1)
			csiAf = -1;
		if (csiAi > +1)
			csiAi = +1;
		if (csiAi < -1)
			csiAi = -1;
		if (csiBf > +1)
			csiBf = +1;
		if (csiBf < -1)
			csiBf = -1;
		if (csiBi > +1)
			csiBi = +1;
		if (csiBi < -1)
			csiBi = -1;
		csi_A = 0.5*(csiAi + csiAf);
		csi_B = 0.5*(csiBi + csiBf);
		if (inverted == true)
			csi_B = -csi_B;
	}
	//Tensor rotação no csiA:
	Matrix Q_Acsi = *surf1->Q_AAic*(0.5 - 0.5*csi_A) + *surf1->Q_BAic*(0.5 + 0.5*csi_A);
	//Tensor rotação no csiB:
	Matrix Q_Bcsi = *surf2->Q_AAic*(0.5 - 0.5*csi_B) + *surf2->Q_BAic*(0.5 + 0.5*csi_B);
	//Posição no csiA:
	Matrix x_Acsi = (*surf1->x_AAi+uAA)*(0.5 - 0.5*csi_A) + (*surf1->x_BAi+uBA)*(0.5 + 0.5*csi_A);
	//Posição no csiB:
	Matrix x_Bcsi = (*surf2->x_AAi+uAB)*(0.5 - 0.5*csi_B) + (*surf2->x_BAi+uBB)*(0.5 + 0.5*csi_B);
	//Vetor e1_A atualizado (em csi)
	Matrix e1_A = Q_Acsi*(*db.CS[surf1->csA - 1]->E1);
	Matrix e2_A = Q_Acsi*(*db.CS[surf1->csA - 1]->E2);
	Matrix e3_A = Q_Acsi*(*db.CS[surf1->csA - 1]->E3);
	//Vetor e1_B atualizado (em csi)
	Matrix e1_B = Q_Bcsi*(*db.CS[surf2->csB - 1]->E1);
	Matrix e2_B = Q_Bcsi*(*db.CS[surf2->csB - 1]->E2);
	Matrix e3_B = Q_Bcsi*(*db.CS[surf2->csB - 1]->E3);
	//Estimativa dos theta_A e theta_B
	double theta_A = 0;
	double theta_B = 0;
	double ndAB = norm(x_Acsi - x_Bcsi);
	Matrix dAB = Matrix(3);
	if (ndAB != 0.0)
		dAB = (1.0 / ndAB)*(x_Acsi - x_Bcsi);
	Matrix dAB_planeA = -1.0*dAB - dot(-1.0*dAB, e3_A)*e3_A;
	Matrix dAB_planeB = dAB - dot(dAB, e3_B)*e3_B;
	double cos_e1_A = dot(dAB_planeA, e1_A);
	double cos_e2_A = dot(dAB_planeA, e2_A);
	//Evitando erros de arredondamento - cos
	if (cos_e1_A > 1.0)
		cos_e1_A = 1.0;
	if (cos_e1_A < -1.0)
		cos_e1_A = -1.0;
	if (cos_e2_A > 1.0)
		cos_e2_A = 1.0;
	if (cos_e2_A < -1.0)
		cos_e2_A = -1.0;
	//Primeiro ou quarto quadrantes
	if (cos_e1_A >= 0.0)
	{
		if (cos_e2_A >= 0.0)//primeiro
			theta_A = abs(acos(cos_e1_A));
		else //quarto
			theta_A = -abs(acos(cos_e1_A));
	}
	//Segundo ou terceiro quadrantes
	else
	{
		if (cos_e2_A >= 0.0)//segundo
			theta_A = abs(acos(cos_e2_A)) + PI / 2;
		else //terceiro
			theta_A = abs(acos(-cos_e1_A)) + PI;
	}
	double cos_e1_B = dot(dAB_planeB, e1_B);
	double cos_e2_B = dot(dAB_planeB, e2_B);
	//Evitando erros de arredondamento - cos
	if (cos_e1_B > 1.0)
		cos_e1_B = 1.0;
	if (cos_e1_B < -1.0)
		cos_e1_B = -1.0;
	if (cos_e2_B > 1.0)
		cos_e2_B = 1.0;
	if (cos_e2_B < -1.0)
		cos_e2_B = -1.0;
	//Primeiro ou quarto quadrantes
	if (cos_e1_B >= 0.0)
	{
		if (cos_e2_B >= 0.0)//primeiro
			theta_B = abs(acos(cos_e1_B));
		else //quarto
			theta_B = -abs(acos(cos_e1_B));
	}
	//Segundo ou terceiro quadrantes
	else
	{
		if (cos_e2_B >= 0.0)//segundo
			theta_B = abs(acos(cos_e2_B)) + PI / 2;
		else //terceiro
			theta_B = abs(acos(-cos_e1_B)) + PI;
	}
	for (int ip = 0; ip < c_data->n_solutions; ip++)
	{
		//Preenchendo as coordenadas convectivas:
		c_data->convective[ip][0] = csi_A;
		c_data->convective[ip][1] = theta_A;
		c_data->convective[ip][2] = csi_B;
		c_data->convective[ip][3] = theta_B;
	}
	
}

//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
double FlexibleSECylinder_1_FlexibleSECylinder_1::ObjectivePhase1(Matrix& mc)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
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
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	aA = &surf1->a;
	aB = &surf2->a;
	bA = &surf1->b;
	bB = &surf2->b;
	eA = &surf1->e;
	eB = &surf2->e;
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_A->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_A->getMatrix();
	dduiB = surf2->ddui_A->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_AAi->getMatrix();
	xBBi = surf2->x_BAi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	QBBi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
		QBBi[i] = new double[3];
	}
	//Salvando variáveis locais para montagem de superfícies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_AAi->MatrixToPtr(QABi, 3);
	surf2->Q_BAi->MatrixToPtr(QBBi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double v[600];
	double *c = mc.getMatrix();
	double Ob;

	v[411] = 2e0 / (*eB);
	v[410] = 0.5e0*(1e0 + c[2]);
	v[409] = 0.5e0 - 0.5e0*c[2];
	v[404] = 2e0 / (*eA);
	v[403] = 0.5e0*(1e0 + c[0]);
	v[402] = 0.5e0 - 0.5e0*c[0];
	v[397] = Power(dB[11], 2);
	v[396] = 0.5e0*dB[11] * dB[9];
	v[395] = 0.5e0*dB[10];
	v[394] = Power(dB[10], 2);
	v[407] = v[394] + v[397];
	v[393] = dB[9] * v[395];
	v[392] = Power(dB[9], 2);
	v[391] = Power(dB[5], 2);
	v[390] = 0.5e0*dB[3] * dB[5];
	v[389] = 0.5e0*dB[4];
	v[388] = Power(dB[4], 2);
	v[405] = v[388] + v[391];
	v[387] = dB[3] * v[389];
	v[386] = Power(dB[3], 2);
	v[385] = Power(dA[11], 2);
	v[384] = 0.5e0*dA[11] * dA[9];
	v[383] = 0.5e0*dA[10];
	v[382] = Power(dA[10], 2);
	v[400] = v[382] + v[385];
	v[381] = dA[9] * v[383];
	v[380] = Power(dA[9], 2);
	v[379] = Power(dA[5], 2);
	v[378] = 0.5e0*dA[3] * dA[5];
	v[377] = 0.5e0*dA[4];
	v[376] = Power(dA[4], 2);
	v[398] = v[376] + v[379];
	v[375] = dA[3] * v[377];
	v[374] = Power(dA[3], 2);
	v[373] = cos(c[3]);
	v[372] = sin(c[3]);
	v[371] = cos(c[1]);
	v[370] = sin(c[1]);
	v[126] = dA[5] * v[377];
	v[142] = dA[11] * v[383];
	v[193] = dB[5] * v[389];
	v[209] = dB[11] * v[395];
	v[113] = 4e0 / (4e0 + v[374] + v[398]);
	v[399] = -0.5e0*v[113];
	v[116] = 1e0 + v[398] * v[399];
	v[117] = v[113] * (-dA[5] + v[375]);
	v[118] = v[113] * (dA[4] + v[378]);
	v[120] = v[113] * (dA[5] + v[375]);
	v[122] = 1e0 + (v[374] + v[379])*v[399];
	v[123] = v[113] * (-dA[3] + v[126]);
	v[125] = v[113] * (-dA[4] + v[378]);
	v[127] = v[113] * (dA[3] + v[126]);
	v[128] = 1e0 + (v[374] + v[376])*v[399];
	v[129] = 4e0 / (4e0 + v[380] + v[400]);
	v[401] = -0.5e0*v[129];
	v[132] = 1e0 + v[400] * v[401];
	v[133] = v[129] * (-dA[11] + v[381]);
	v[134] = v[129] * (dA[10] + v[384]);
	v[136] = v[129] * (dA[11] + v[381]);
	v[138] = 1e0 + (v[380] + v[385])*v[401];
	v[139] = v[129] * (-dA[9] + v[142]);
	v[141] = v[129] * (-dA[10] + v[384]);
	v[143] = v[129] * (dA[9] + v[142]);
	v[144] = 1e0 + (v[380] + v[382])*v[401];
	v[172] = 1e0 / Power(Power(fabs(v[370]), v[404]) + Power(fabs(v[371]), v[404]), 1 / v[404]);
	v[174] = (*aA)*v[172] * v[371];
	v[176] = (*bA)*v[172] * v[370];
	v[180] = 4e0 / (4e0 + v[386] + v[405]);
	v[406] = -0.5e0*v[180];
	v[183] = 1e0 + v[405] * v[406];
	v[184] = v[180] * (-dB[5] + v[387]);
	v[185] = v[180] * (dB[4] + v[390]);
	v[187] = v[180] * (dB[5] + v[387]);
	v[189] = 1e0 + (v[386] + v[391])*v[406];
	v[190] = v[180] * (-dB[3] + v[193]);
	v[192] = v[180] * (-dB[4] + v[390]);
	v[194] = v[180] * (dB[3] + v[193]);
	v[195] = 1e0 + (v[386] + v[388])*v[406];
	v[196] = 4e0 / (4e0 + v[392] + v[407]);
	v[408] = -0.5e0*v[196];
	v[199] = 1e0 + v[407] * v[408];
	v[200] = v[196] * (-dB[11] + v[393]);
	v[201] = v[196] * (dB[10] + v[396]);
	v[203] = v[196] * (dB[11] + v[393]);
	v[205] = 1e0 + (v[392] + v[397])*v[408];
	v[206] = v[196] * (-dB[9] + v[209]);
	v[208] = v[196] * (-dB[10] + v[396]);
	v[210] = v[196] * (dB[9] + v[209]);
	v[211] = 1e0 + (v[392] + v[394])*v[408];
	v[239] = 1e0 / Power(Power(fabs(v[372]), v[411]) + Power(fabs(v[373]), v[411]), 1 / v[411]);
	v[241] = (*aB)*v[239] * v[373];
	v[243] = (*bB)*v[239] * v[372];
	(Ob) = 0.5e0*(Power(v[402] * (dA[0] + (QAAi[0][0] * v[116] + QAAi[1][0] * v[117] + QAAi[2][0] * v[118])*v[174] +
		(QAAi[0][1] * v[116] + QAAi[1][1] * v[117] + QAAi[2][1] * v[118])*v[176] + xAAi[0]) - v[409] * (dB[0] +
		(QABi[0][0] * v[183] + QABi[1][0] * v[184] + QABi[2][0] * v[185])*v[241] + (QABi[0][1] * v[183] + QABi[1][1] * v[184]
		+ QABi[2][1] * v[185])*v[243] + xABi[0]) + v[403] * (dA[6] + (QBAi[0][0] * v[132] + QBAi[1][0] * v[133]
		+ QBAi[2][0] * v[134])*v[174] + (QBAi[0][1] * v[132] + QBAi[1][1] * v[133] + QBAi[2][1] * v[134])*v[176] + xBAi[0])
		- v[410] * (dB[6] + (QBBi[0][0] * v[199] + QBBi[1][0] * v[200] + QBBi[2][0] * v[201])*v[241] + (QBBi[0][1] * v[199]
		+ QBBi[1][1] * v[200] + QBBi[2][1] * v[201])*v[243] + xBBi[0]), 2) + Power(v[402] * (dA[1] + (QAAi[0][0] * v[120]
		+ QAAi[1][0] * v[122] + QAAi[2][0] * v[123])*v[174] + (QAAi[0][1] * v[120] + QAAi[1][1] * v[122] + QAAi[2][1] * v[123]
		)*v[176] + xAAi[1]) - v[409] * (dB[1] + (QABi[0][0] * v[187] + QABi[1][0] * v[189] + QABi[2][0] * v[190])*v[241] +
		(QABi[0][1] * v[187] + QABi[1][1] * v[189] + QABi[2][1] * v[190])*v[243] + xABi[1]) + v[403] * (dA[7] +
		(QBAi[0][0] * v[136] + QBAi[1][0] * v[138] + QBAi[2][0] * v[139])*v[174] + (QBAi[0][1] * v[136] + QBAi[1][1] * v[138]
		+ QBAi[2][1] * v[139])*v[176] + xBAi[1]) - v[410] * (dB[7] + (QBBi[0][0] * v[203] + QBBi[1][0] * v[205]
		+ QBBi[2][0] * v[206])*v[241] + (QBBi[0][1] * v[203] + QBBi[1][1] * v[205] + QBBi[2][1] * v[206])*v[243] + xBBi[1])
		, 2) + Power(v[402] * (dA[2] + (QAAi[0][0] * v[125] + QAAi[1][0] * v[127] + QAAi[2][0] * v[128])*v[174] +
		(QAAi[0][1] * v[125] + QAAi[1][1] * v[127] + QAAi[2][1] * v[128])*v[176] + xAAi[2]) - v[409] * (dB[2] +
		(QABi[0][0] * v[192] + QABi[1][0] * v[194] + QABi[2][0] * v[195])*v[241] + (QABi[0][1] * v[192] + QABi[1][1] * v[194]
		+ QABi[2][1] * v[195])*v[243] + xABi[2]) + v[403] * (dA[8] + (QBAi[0][0] * v[141] + QBAi[1][0] * v[143]
		+ QBAi[2][0] * v[144])*v[174] + (QBAi[0][1] * v[141] + QBAi[1][1] * v[143] + QBAi[2][1] * v[144])*v[176] + xBAi[2])
		- v[410] * (dB[8] + (QBBi[0][0] * v[208] + QBBi[1][0] * v[210] + QBBi[2][0] * v[211])*v[241] + (QBBi[0][1] * v[208]
		+ QBBi[1][1] * v[210] + QBBi[2][1] * v[211])*v[243] + xBBi[2]), 2));
	
	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
		delete[] QBBi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
	delete[] QBBi;

	return Ob;
}

//Calcula o Gradiente da função objetivo - Phase 1
void FlexibleSECylinder_1_FlexibleSECylinder_1::GradientPhase1(Matrix& mc, Matrix& mGra)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
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
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	aA = &surf1->a;
	aB = &surf2->a;
	bA = &surf1->b;
	bB = &surf2->b;
	eA = &surf1->e;
	eB = &surf2->e;
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_A->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_A->getMatrix();
	dduiB = surf2->ddui_A->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_AAi->getMatrix();
	xBBi = surf2->x_BAi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	QBBi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
		QBBi[i] = new double[3];
	}
	//Salvando variáveis locais para montagem de superfícies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_AAi->MatrixToPtr(QABi, 3);
	surf2->Q_BAi->MatrixToPtr(QBBi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double v[600];
	double *c = mc.getMatrix();
	double Gra[4];

	v[415] = 2e0 / (*eB);
	v[417] = 1e0 / v[415];
	v[416] = -1e0 + v[415];
	v[414] = 0.5e0*(1e0 + c[2]);
	v[413] = 0.5e0 - 0.5e0*c[2];
	v[404] = 2e0 / (*eA);
	v[406] = 1e0 / v[404];
	v[405] = -1e0 + v[404];
	v[403] = 0.5e0*(1e0 + c[0]);
	v[402] = 0.5e0 - 0.5e0*c[0];
	v[397] = Power(dB[11], 2);
	v[396] = 0.5e0*dB[11] * dB[9];
	v[395] = 0.5e0*dB[10];
	v[394] = Power(dB[10], 2);
	v[411] = v[394] + v[397];
	v[393] = dB[9] * v[395];
	v[392] = Power(dB[9], 2);
	v[391] = Power(dB[5], 2);
	v[390] = 0.5e0*dB[3] * dB[5];
	v[389] = 0.5e0*dB[4];
	v[388] = Power(dB[4], 2);
	v[409] = v[388] + v[391];
	v[387] = dB[3] * v[389];
	v[386] = Power(dB[3], 2);
	v[385] = Power(dA[11], 2);
	v[384] = 0.5e0*dA[11] * dA[9];
	v[383] = 0.5e0*dA[10];
	v[382] = Power(dA[10], 2);
	v[400] = v[382] + v[385];
	v[381] = dA[9] * v[383];
	v[380] = Power(dA[9], 2);
	v[379] = Power(dA[5], 2);
	v[378] = 0.5e0*dA[3] * dA[5];
	v[377] = 0.5e0*dA[4];
	v[376] = Power(dA[4], 2);
	v[398] = v[376] + v[379];
	v[375] = dA[3] * v[377];
	v[374] = Power(dA[3], 2);
	v[373] = cos(c[3]);
	v[372] = sin(c[3]);
	v[371] = cos(c[1]);
	v[370] = sin(c[1]);
	v[253] = fabs(v[370]);
	v[251] = fabs(v[371]);
	v[260] = fabs(v[372]);
	v[258] = fabs(v[373]);
	v[126] = dA[5] * v[377];
	v[142] = dA[11] * v[383];
	v[193] = dB[5] * v[389];
	v[209] = dB[11] * v[395];
	v[113] = 4e0 / (4e0 + v[374] + v[398]);
	v[399] = -0.5e0*v[113];
	v[116] = 1e0 + v[398] * v[399];
	v[117] = v[113] * (-dA[5] + v[375]);
	v[118] = v[113] * (dA[4] + v[378]);
	v[120] = v[113] * (dA[5] + v[375]);
	v[122] = 1e0 + (v[374] + v[379])*v[399];
	v[123] = v[113] * (-dA[3] + v[126]);
	v[125] = v[113] * (-dA[4] + v[378]);
	v[127] = v[113] * (dA[3] + v[126]);
	v[128] = 1e0 + (v[374] + v[376])*v[399];
	v[129] = 4e0 / (4e0 + v[380] + v[400]);
	v[401] = -0.5e0*v[129];
	v[132] = 1e0 + v[400] * v[401];
	v[133] = v[129] * (-dA[11] + v[381]);
	v[134] = v[129] * (dA[10] + v[384]);
	v[136] = v[129] * (dA[11] + v[381]);
	v[138] = 1e0 + (v[380] + v[385])*v[401];
	v[139] = v[129] * (-dA[9] + v[142]);
	v[141] = v[129] * (-dA[10] + v[384]);
	v[143] = v[129] * (dA[9] + v[142]);
	v[144] = 1e0 + (v[380] + v[382])*v[401];
	v[145] = QAAi[0][0] * v[116] + QAAi[1][0] * v[117] + QAAi[2][0] * v[118];
	v[146] = QAAi[0][1] * v[116] + QAAi[1][1] * v[117] + QAAi[2][1] * v[118];
	v[148] = QAAi[0][0] * v[120] + QAAi[1][0] * v[122] + QAAi[2][0] * v[123];
	v[149] = QAAi[0][1] * v[120] + QAAi[1][1] * v[122] + QAAi[2][1] * v[123];
	v[151] = QAAi[0][0] * v[125] + QAAi[1][0] * v[127] + QAAi[2][0] * v[128];
	v[152] = QAAi[0][1] * v[125] + QAAi[1][1] * v[127] + QAAi[2][1] * v[128];
	v[154] = QBAi[0][0] * v[132] + QBAi[1][0] * v[133] + QBAi[2][0] * v[134];
	v[155] = QBAi[0][1] * v[132] + QBAi[1][1] * v[133] + QBAi[2][1] * v[134];
	v[157] = QBAi[0][0] * v[136] + QBAi[1][0] * v[138] + QBAi[2][0] * v[139];
	v[158] = QBAi[0][1] * v[136] + QBAi[1][1] * v[138] + QBAi[2][1] * v[139];
	v[160] = QBAi[0][0] * v[141] + QBAi[1][0] * v[143] + QBAi[2][0] * v[144];
	v[161] = QBAi[0][1] * v[141] + QBAi[1][1] * v[143] + QBAi[2][1] * v[144];
	v[250] = Power(v[251], v[404]) + Power(v[253], v[404]);
	v[254] = Power(v[250], -1e0 - v[406])*(-(Power(v[253], v[405])*v[371] * _copysign(1.e0, v[370])) + Power
		(v[251], v[405])*v[370] * _copysign(1.e0, v[371]));
	v[172] = 1e0 / Power(v[250], v[406]);
	v[408] = v[172] * v[370];
	v[407] = v[172] * v[371];
	v[256] = (*bA)*(v[254] * v[370] + v[407]);
	v[255] = (*aA)*(v[254] * v[371] - v[408]);
	v[174] = (*aA)*v[407];
	v[176] = (*bA)*v[408];
	v[271] = dA[6] + v[154] * v[174] + v[155] * v[176] + xBAi[0];
	v[270] = dA[0] + v[145] * v[174] + v[146] * v[176] + xAAi[0];
	v[268] = dA[7] + v[157] * v[174] + v[158] * v[176] + xBAi[1];
	v[267] = dA[1] + v[148] * v[174] + v[149] * v[176] + xAAi[1];
	v[265] = dA[8] + v[160] * v[174] + v[161] * v[176] + xBAi[2];
	v[264] = dA[2] + v[151] * v[174] + v[152] * v[176] + xAAi[2];
	v[180] = 4e0 / (4e0 + v[386] + v[409]);
	v[410] = -0.5e0*v[180];
	v[183] = 1e0 + v[409] * v[410];
	v[184] = v[180] * (-dB[5] + v[387]);
	v[185] = v[180] * (dB[4] + v[390]);
	v[187] = v[180] * (dB[5] + v[387]);
	v[189] = 1e0 + (v[386] + v[391])*v[410];
	v[190] = v[180] * (-dB[3] + v[193]);
	v[192] = v[180] * (-dB[4] + v[390]);
	v[194] = v[180] * (dB[3] + v[193]);
	v[195] = 1e0 + (v[386] + v[388])*v[410];
	v[196] = 4e0 / (4e0 + v[392] + v[411]);
	v[412] = -0.5e0*v[196];
	v[199] = 1e0 + v[411] * v[412];
	v[200] = v[196] * (-dB[11] + v[393]);
	v[201] = v[196] * (dB[10] + v[396]);
	v[203] = v[196] * (dB[11] + v[393]);
	v[205] = 1e0 + (v[392] + v[397])*v[412];
	v[206] = v[196] * (-dB[9] + v[209]);
	v[208] = v[196] * (-dB[10] + v[396]);
	v[210] = v[196] * (dB[9] + v[209]);
	v[211] = 1e0 + (v[392] + v[394])*v[412];
	v[212] = QABi[0][0] * v[183] + QABi[1][0] * v[184] + QABi[2][0] * v[185];
	v[213] = QABi[0][1] * v[183] + QABi[1][1] * v[184] + QABi[2][1] * v[185];
	v[215] = QABi[0][0] * v[187] + QABi[1][0] * v[189] + QABi[2][0] * v[190];
	v[216] = QABi[0][1] * v[187] + QABi[1][1] * v[189] + QABi[2][1] * v[190];
	v[218] = QABi[0][0] * v[192] + QABi[1][0] * v[194] + QABi[2][0] * v[195];
	v[219] = QABi[0][1] * v[192] + QABi[1][1] * v[194] + QABi[2][1] * v[195];
	v[221] = QBBi[0][0] * v[199] + QBBi[1][0] * v[200] + QBBi[2][0] * v[201];
	v[222] = QBBi[0][1] * v[199] + QBBi[1][1] * v[200] + QBBi[2][1] * v[201];
	v[224] = QBBi[0][0] * v[203] + QBBi[1][0] * v[205] + QBBi[2][0] * v[206];
	v[225] = QBBi[0][1] * v[203] + QBBi[1][1] * v[205] + QBBi[2][1] * v[206];
	v[227] = QBBi[0][0] * v[208] + QBBi[1][0] * v[210] + QBBi[2][0] * v[211];
	v[228] = QBBi[0][1] * v[208] + QBBi[1][1] * v[210] + QBBi[2][1] * v[211];
	v[257] = Power(v[258], v[415]) + Power(v[260], v[415]);
	v[261] = Power(v[257], -1e0 - v[417])*(-(Power(v[260], v[416])*v[373] * _copysign(1.e0, v[372])) + Power
		(v[258], v[416])*v[372] * _copysign(1.e0, v[373]));
	v[239] = 1e0 / Power(v[257], v[417]);
	v[419] = v[239] * v[372];
	v[418] = v[239] * v[373];
	v[263] = (*bB)*(v[261] * v[372] + v[418]);
	v[262] = (*aB)*(v[261] * v[373] - v[419]);
	v[241] = (*aB)*v[418];
	v[243] = (*bB)*v[419];
	v[283] = dB[6] + v[221] * v[241] + v[222] * v[243] + xBBi[0];
	v[282] = dB[0] + v[212] * v[241] + v[213] * v[243] + xABi[0];
	v[280] = dB[7] + v[224] * v[241] + v[225] * v[243] + xBBi[1];
	v[279] = dB[1] + v[215] * v[241] + v[216] * v[243] + xABi[1];
	v[277] = dB[8] + v[227] * v[241] + v[228] * v[243] + xBBi[2];
	v[276] = dB[2] + v[218] * v[241] + v[219] * v[243] + xABi[2];
	v[290] = v[270] * v[402] + v[271] * v[403] - v[282] * v[413] - v[283] * v[414];
	v[422] = 2e0*v[290];
	v[289] = v[267] * v[402] + v[268] * v[403] - v[279] * v[413] - v[280] * v[414];
	v[421] = 2e0*v[289];
	v[288] = v[264] * v[402] + v[265] * v[403] - v[276] * v[413] - v[277] * v[414];
	v[420] = 2e0*v[288];
	Gra[0] = 0.5e0*((-0.5e0*v[264] + 0.5e0*v[265])*v[420] + (-0.5e0*v[267] + 0.5e0*v[268])*v[421] + (
		-0.5e0*v[270] + 0.5e0*v[271])*v[422]);
	Gra[1] = 0.5e0*(((v[151] * v[255] + v[152] * v[256])*v[402] + (v[160] * v[255] + v[161] * v[256])*v[403])*v[420] + (
		(v[148] * v[255] + v[149] * v[256])*v[402] + (v[157] * v[255] + v[158] * v[256])*v[403])*v[421] + ((v[145] * v[255]
		+ v[146] * v[256])*v[402] + (v[154] * v[255] + v[155] * v[256])*v[403])*v[422]);
	Gra[2] = 0.5e0*((0.5e0*v[276] - 0.5e0*v[277])*v[420] + (0.5e0*v[279] - 0.5e0*v[280])*v[421] + (0.5e0*v[282]
		- 0.5e0*v[283])*v[422]);
	Gra[3] = 0.5e0*((-((v[218] * v[262] + v[219] * v[263])*v[413]) - (v[227] * v[262] + v[228] * v[263])*v[414])*v[420]
		+ (-((v[215] * v[262] + v[216] * v[263])*v[413]) - (v[224] * v[262] + v[225] * v[263])*v[414])*v[421] + (-(
		(v[212] * v[262] + v[213] * v[263])*v[413]) - (v[221] * v[262] + v[222] * v[263])*v[414])*v[422]);

	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
		delete[] QBBi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
	delete[] QBBi;
}

//Calcula a Hessiana da função objetivo - Phase 1
void FlexibleSECylinder_1_FlexibleSECylinder_1::HessianPhase1(Matrix& mc, Matrix& mHes)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
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
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	aA = &surf1->a;
	aB = &surf2->a;
	bA = &surf1->b;
	bB = &surf2->b;
	eA = &surf1->e;
	eB = &surf2->e;
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_A->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_A->getMatrix();
	dduiB = surf2->ddui_A->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_AAi->getMatrix();
	xBBi = surf2->x_BAi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	QBBi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
		QBBi[i] = new double[3];
	}
	//Salvando variáveis locais para montagem de superfícies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_AAi->MatrixToPtr(QABi, 3);
	surf2->Q_BAi->MatrixToPtr(QBBi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double v[450];
	double *c = mc.getMatrix();
	double Hes[4][4];
	//Pontos de singularidade - função phi
	double test1 = 2 * c[1] / PI;
	test1 = abs(test1 - int(test1));//testa o quão próximo de inteiro é esse número
	if (test1 < tol_small_1)
		c[1] += tol_small_1;
	double test3 = 2 * c[3] / PI;
	test3 = abs(test3 - int(test3));//testa o quão próximo de inteiro é esse número
	if (test3 < tol_small_1)
		c[3] += tol_small_1;
	

	int i01; int i02;
	v[417] = 2e0 / (*eB);
	v[419] = 1e0 / v[417];
	v[418] = -1e0 + v[417];
	v[416] = 0.5e0*(1e0 + c[2]);
	v[415] = 0.5e0 - 0.5e0*c[2];
	v[405] = 2e0 / (*eA);
	v[407] = 1e0 / v[405];
	v[406] = -1e0 + v[405];
	v[404] = 0.5e0*(1e0 + c[0]);
	v[403] = 0.5e0 - 0.5e0*c[0];
	v[398] = Power(dB[11], 2);
	v[397] = 0.5e0*dB[11] * dB[9];
	v[396] = 0.5e0*dB[10];
	v[395] = Power(dB[10], 2);
	v[413] = v[395] + v[398];
	v[394] = dB[9] * v[396];
	v[393] = Power(dB[9], 2);
	v[392] = Power(dB[5], 2);
	v[391] = 0.5e0*dB[3] * dB[5];
	v[390] = 0.5e0*dB[4];
	v[389] = Power(dB[4], 2);
	v[411] = v[389] + v[392];
	v[388] = dB[3] * v[390];
	v[387] = Power(dB[3], 2);
	v[386] = Power(dA[11], 2);
	v[385] = 0.5e0*dA[11] * dA[9];
	v[384] = 0.5e0*dA[10];
	v[383] = Power(dA[10], 2);
	v[401] = v[383] + v[386];
	v[382] = dA[9] * v[384];
	v[381] = Power(dA[9], 2);
	v[380] = Power(dA[5], 2);
	v[379] = 0.5e0*dA[3] * dA[5];
	v[378] = 0.5e0*dA[4];
	v[377] = Power(dA[4], 2);
	v[399] = v[377] + v[380];
	v[376] = dA[3] * v[378];
	v[375] = Power(dA[3], 2);
	v[374] = cos(c[3]);
	v[373] = sin(c[3]);
	v[372] = cos(c[1]);
	v[371] = sin(c[1]);
	v[296] = _copysign(1.e0, v[371]);
	v[297] = v[296] * v[372];
	v[253] = fabs(v[371]);
	v[298] = _copysign(1.e0, v[372]);
	v[299] = -(v[298] * v[371]);
	v[251] = fabs(v[372]);
	v[301] = _copysign(1.e0, v[373]);
	v[302] = v[301] * v[374];
	v[260] = fabs(v[373]);
	v[303] = _copysign(1.e0, v[374]);
	v[304] = -(v[303] * v[373]);
	v[258] = fabs(v[374]);
	v[126] = dA[5] * v[378];
	v[142] = dA[11] * v[384];
	v[193] = dB[5] * v[390];
	v[209] = dB[11] * v[396];
	v[113] = 4e0 / (4e0 + v[375] + v[399]);
	v[400] = -0.5e0*v[113];
	v[116] = 1e0 + v[399] * v[400];
	v[117] = v[113] * (-dA[5] + v[376]);
	v[118] = v[113] * (dA[4] + v[379]);
	v[120] = v[113] * (dA[5] + v[376]);
	v[122] = 1e0 + (v[375] + v[380])*v[400];
	v[123] = v[113] * (-dA[3] + v[126]);
	v[125] = v[113] * (-dA[4] + v[379]);
	v[127] = v[113] * (dA[3] + v[126]);
	v[128] = 1e0 + (v[375] + v[377])*v[400];
	v[129] = 4e0 / (4e0 + v[381] + v[401]);
	v[402] = -0.5e0*v[129];
	v[132] = 1e0 + v[401] * v[402];
	v[133] = v[129] * (-dA[11] + v[382]);
	v[134] = v[129] * (dA[10] + v[385]);
	v[136] = v[129] * (dA[11] + v[382]);
	v[138] = 1e0 + (v[381] + v[386])*v[402];
	v[139] = v[129] * (-dA[9] + v[142]);
	v[141] = v[129] * (-dA[10] + v[385]);
	v[143] = v[129] * (dA[9] + v[142]);
	v[144] = 1e0 + (v[381] + v[383])*v[402];
	v[145] = QAAi[0][0] * v[116] + QAAi[1][0] * v[117] + QAAi[2][0] * v[118];
	v[146] = QAAi[0][1] * v[116] + QAAi[1][1] * v[117] + QAAi[2][1] * v[118];
	v[148] = QAAi[0][0] * v[120] + QAAi[1][0] * v[122] + QAAi[2][0] * v[123];
	v[149] = QAAi[0][1] * v[120] + QAAi[1][1] * v[122] + QAAi[2][1] * v[123];
	v[151] = QAAi[0][0] * v[125] + QAAi[1][0] * v[127] + QAAi[2][0] * v[128];
	v[152] = QAAi[0][1] * v[125] + QAAi[1][1] * v[127] + QAAi[2][1] * v[128];
	v[154] = QBAi[0][0] * v[132] + QBAi[1][0] * v[133] + QBAi[2][0] * v[134];
	v[155] = QBAi[0][1] * v[132] + QBAi[1][1] * v[133] + QBAi[2][1] * v[134];
	v[157] = QBAi[0][0] * v[136] + QBAi[1][0] * v[138] + QBAi[2][0] * v[139];
	v[158] = QBAi[0][1] * v[136] + QBAi[1][1] * v[138] + QBAi[2][1] * v[139];
	v[160] = QBAi[0][0] * v[141] + QBAi[1][0] * v[143] + QBAi[2][0] * v[144];
	v[161] = QBAi[0][1] * v[141] + QBAi[1][1] * v[143] + QBAi[2][1] * v[144];
	v[312] = -1e0 - 1e0 / v[405];
	v[310] = Power(v[253], v[406]);
	v[308] = Power(v[251], v[406]);
	v[305] = (v[299] * v[308] + v[297] * v[310])*v[405];
	v[311] = -1e0 + v[406];
	v[250] = Power(v[251], v[405]) + Power(v[253], v[405]);
	v[306] = Power(v[250], v[312]);
	v[314] = -((Power(v[250], -2e0 - v[407])*(v[305] * v[305])*v[312]) / v[405]) - v[306] * (v[298] * (-(v[308] * v[372]
		) - Power(v[251], v[311])*v[299] * v[371] * v[406]) + v[296] * (-(v[310] * v[371]) + Power(v[253], v[311]
		)*v[297] * v[372] * v[406]));
	v[254] = -((v[305] * v[306]) / v[405]);
	v[408] = 2e0*v[254];
	v[172] = 1e0 / Power(v[250], v[407]);
	v[410] = v[172] * v[371];
	v[409] = v[172] * v[372];
	v[317] = (*aA)*v[409];
	v[319] = -v[317] + (*aA)*(v[314] * v[372] - v[371] * v[408]);
	v[315] = (*bA)*v[410];
	v[316] = -v[315] + (*bA)*(v[314] * v[371] + v[372] * v[408]);
	v[256] = (*bA)*(v[254] * v[371] + v[409]);
	v[255] = (*aA)*(v[254] * v[372] - v[410]);
	v[325] = v[151] * v[255] + v[152] * v[256];
	v[324] = v[160] * v[255] + v[161] * v[256];
	v[323] = v[148] * v[255] + v[149] * v[256];
	v[322] = v[157] * v[255] + v[158] * v[256];
	v[321] = v[145] * v[255] + v[146] * v[256];
	v[320] = v[154] * v[255] + v[155] * v[256];
	v[275] = v[321] * v[403] + v[320] * v[404];
	v[425] = 2e0*v[275];
	v[274] = v[323] * v[403] + v[322] * v[404];
	v[424] = 2e0*v[274];
	v[273] = v[325] * v[403] + v[324] * v[404];
	v[423] = 2e0*v[273];
	v[271] = dA[6] + v[155] * v[315] + v[154] * v[317] + xBAi[0];
	v[270] = dA[0] + v[146] * v[315] + v[145] * v[317] + xAAi[0];
	v[272] = -0.5e0*v[270] + 0.5e0*v[271];
	v[268] = dA[7] + v[158] * v[315] + v[157] * v[317] + xBAi[1];
	v[267] = dA[1] + v[149] * v[315] + v[148] * v[317] + xAAi[1];
	v[269] = -0.5e0*v[267] + 0.5e0*v[268];
	v[265] = dA[8] + v[161] * v[315] + v[160] * v[317] + xBAi[2];
	v[264] = dA[2] + v[152] * v[315] + v[151] * v[317] + xAAi[2];
	v[266] = -0.5e0*v[264] + 0.5e0*v[265];
	v[180] = 4e0 / (4e0 + v[387] + v[411]);
	v[412] = -0.5e0*v[180];
	v[183] = 1e0 + v[411] * v[412];
	v[184] = v[180] * (-dB[5] + v[388]);
	v[185] = v[180] * (dB[4] + v[391]);
	v[187] = v[180] * (dB[5] + v[388]);
	v[189] = 1e0 + (v[387] + v[392])*v[412];
	v[190] = v[180] * (-dB[3] + v[193]);
	v[192] = v[180] * (-dB[4] + v[391]);
	v[194] = v[180] * (dB[3] + v[193]);
	v[195] = 1e0 + (v[387] + v[389])*v[412];
	v[196] = 4e0 / (4e0 + v[393] + v[413]);
	v[414] = -0.5e0*v[196];
	v[199] = 1e0 + v[413] * v[414];
	v[200] = v[196] * (-dB[11] + v[394]);
	v[201] = v[196] * (dB[10] + v[397]);
	v[203] = v[196] * (dB[11] + v[394]);
	v[205] = 1e0 + (v[393] + v[398])*v[414];
	v[206] = v[196] * (-dB[9] + v[209]);
	v[208] = v[196] * (-dB[10] + v[397]);
	v[210] = v[196] * (dB[9] + v[209]);
	v[211] = 1e0 + (v[393] + v[395])*v[414];
	v[212] = QABi[0][0] * v[183] + QABi[1][0] * v[184] + QABi[2][0] * v[185];
	v[213] = QABi[0][1] * v[183] + QABi[1][1] * v[184] + QABi[2][1] * v[185];
	v[215] = QABi[0][0] * v[187] + QABi[1][0] * v[189] + QABi[2][0] * v[190];
	v[216] = QABi[0][1] * v[187] + QABi[1][1] * v[189] + QABi[2][1] * v[190];
	v[218] = QABi[0][0] * v[192] + QABi[1][0] * v[194] + QABi[2][0] * v[195];
	v[219] = QABi[0][1] * v[192] + QABi[1][1] * v[194] + QABi[2][1] * v[195];
	v[221] = QBBi[0][0] * v[199] + QBBi[1][0] * v[200] + QBBi[2][0] * v[201];
	v[222] = QBBi[0][1] * v[199] + QBBi[1][1] * v[200] + QBBi[2][1] * v[201];
	v[224] = QBBi[0][0] * v[203] + QBBi[1][0] * v[205] + QBBi[2][0] * v[206];
	v[225] = QBBi[0][1] * v[203] + QBBi[1][1] * v[205] + QBBi[2][1] * v[206];
	v[227] = QBBi[0][0] * v[208] + QBBi[1][0] * v[210] + QBBi[2][0] * v[211];
	v[228] = QBBi[0][1] * v[208] + QBBi[1][1] * v[210] + QBBi[2][1] * v[211];
	v[333] = -1e0 - 1e0 / v[417];
	v[331] = Power(v[260], v[418]);
	v[329] = Power(v[258], v[418]);
	v[326] = (v[304] * v[329] + v[302] * v[331])*v[417];
	v[332] = -1e0 + v[418];
	v[257] = Power(v[258], v[417]) + Power(v[260], v[417]);
	v[327] = Power(v[257], v[333]);
	v[335] = -((Power(v[257], -2e0 - v[419])*(v[326] * v[326])*v[333]) / v[417]) - v[327] * (v[303] * (-(v[329] * v[374]
		) - Power(v[258], v[332])*v[304] * v[373] * v[418]) + v[301] * (-(v[331] * v[373]) + Power(v[260], v[332]
		)*v[302] * v[374] * v[418]));
	v[261] = -((v[326] * v[327]) / v[417]);
	v[420] = 2e0*v[261];
	v[239] = 1e0 / Power(v[257], v[419]);
	v[422] = v[239] * v[373];
	v[421] = v[239] * v[374];
	v[338] = (*aB)*v[421];
	v[340] = -v[338] + (*aB)*(v[335] * v[374] - v[373] * v[420]);
	v[336] = (*bB)*v[422];
	v[337] = -v[336] + (*bB)*(v[335] * v[373] + v[374] * v[420]);
	v[263] = (*bB)*(v[261] * v[373] + v[421]);
	v[262] = (*aB)*(v[261] * v[374] - v[422]);
	v[346] = v[218] * v[262] + v[219] * v[263];
	v[345] = v[227] * v[262] + v[228] * v[263];
	v[344] = v[215] * v[262] + v[216] * v[263];
	v[343] = v[224] * v[262] + v[225] * v[263];
	v[342] = v[212] * v[262] + v[213] * v[263];
	v[341] = v[221] * v[262] + v[222] * v[263];
	v[287] = v[342] * v[415] + v[341] * v[416];
	v[286] = v[344] * v[415] + v[343] * v[416];
	v[285] = v[346] * v[415] + v[345] * v[416];
	v[283] = dB[6] + v[222] * v[336] + v[221] * v[338] + xBBi[0];
	v[282] = dB[0] + v[213] * v[336] + v[212] * v[338] + xABi[0];
	v[284] = -0.5e0*v[282] + 0.5e0*v[283];
	v[431] = 2e0*v[284];
	v[280] = dB[7] + v[225] * v[336] + v[224] * v[338] + xBBi[1];
	v[279] = dB[1] + v[216] * v[336] + v[215] * v[338] + xABi[1];
	v[281] = -0.5e0*v[279] + 0.5e0*v[280];
	v[430] = 2e0*v[281];
	v[277] = dB[8] + v[228] * v[336] + v[227] * v[338] + xBBi[2];
	v[276] = dB[2] + v[219] * v[336] + v[218] * v[338] + xABi[2];
	v[278] = -0.5e0*v[276] + 0.5e0*v[277];
	v[429] = 2e0*v[278];
	v[290] = v[270] * v[403] + v[271] * v[404] - v[282] * v[415] - v[283] * v[416];
	v[426] = 2e0*v[290];
	v[289] = v[267] * v[403] + v[268] * v[404] - v[279] * v[415] - v[280] * v[416];
	v[427] = 2e0*v[289];
	v[288] = v[264] * v[403] + v[265] * v[404] - v[276] * v[415] - v[277] * v[416];
	v[428] = 2e0*v[288];
	Hes[0][0] = 1e0*((v[266] * v[266]) + (v[269] * v[269]) + (v[272] * v[272]));
	Hes[0][1] = 0.5e0*(v[266] * v[423] + v[269] * v[424] + v[272] * v[425] + (0.5e0*v[320] - 0.5e0*v[321])*v[426] +
		(0.5e0*v[322] - 0.5e0*v[323])*v[427] + (0.5e0*v[324] - 0.5e0*v[325])*v[428]);
	Hes[0][2] = 0.5e0*(-(v[266] * v[429]) - v[269] * v[430] - v[272] * v[431]);
	Hes[0][3] = -1e0*(v[266] * v[285] + v[269] * v[286] + v[272] * v[287]);
	Hes[1][1] = 0.5e0*(2e0*(v[273] * v[273]) + 2e0*(v[274] * v[274]) + 2e0*(v[275] * v[275]) + ((v[146] * v[316]
		+ v[145] * v[319])*v[403] + (v[155] * v[316] + v[154] * v[319])*v[404])*v[426] + ((v[149] * v[316] + v[148] * v[319]
		)*v[403] + (v[158] * v[316] + v[157] * v[319])*v[404])*v[427] + ((v[152] * v[316] + v[151] * v[319])*v[403] +
		(v[161] * v[316] + v[160] * v[319])*v[404])*v[428]);
	Hes[1][2] = 0.5e0*(-(v[273] * v[429]) - v[274] * v[430] - v[275] * v[431]);
	Hes[1][3] = 0.5e0*(-(v[285] * v[423]) - v[286] * v[424] - v[287] * v[425]);
	Hes[2][2] = 1e0*((v[278] * v[278]) + (v[281] * v[281]) + (v[284] * v[284]));
	Hes[2][3] = 0.5e0*((-0.5e0*v[341] + 0.5e0*v[342])*v[426] + (-0.5e0*v[343] + 0.5e0*v[344])*v[427] + (
		-0.5e0*v[345] + 0.5e0*v[346])*v[428] + v[285] * v[429] + v[286] * v[430] + v[287] * v[431]);
	Hes[3][3] = 0.5e0*(2e0*(v[285] * v[285]) + 2e0*(v[286] * v[286]) + 2e0*(v[287] * v[287]) + (-((v[213] * v[337]
		+ v[212] * v[340])*v[415]) - (v[222] * v[337] + v[221] * v[340])*v[416])*v[426] + (-((v[216] * v[337]
		+ v[215] * v[340])*v[415]) - (v[225] * v[337] + v[224] * v[340])*v[416])*v[427] + (-((v[219] * v[337]
		+ v[218] * v[340])*v[415]) - (v[228] * v[337] + v[227] * v[340])*v[416])*v[428]);
	for (i01 = 1; i01<4; i01++){
		for (i02 = 0; i02<i01; i02++){
			Hes[i01][i02] = Hes[i02][i01];
		}
	};

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
		mHes(i, j) = Hes[i][j];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
		delete[] QBBi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
	delete[] QBBi;
}

void FlexibleSECylinder_1_FlexibleSECylinder_1::ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft)
{
	double v[30000];		//variável temporária - AceGen
	//Zerando matrizes e vetores
	for (int i = 0; i < 24; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < 24; j++)
			Kc[i][j] = 0.0;
	}
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	FlexibleSECylinder_1* surf1;		//Ponteiro para a superfície 1
	FlexibleSECylinder_1* surf2;		//Ponteiro para a superfície 2
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
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
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;
	surf1 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<FlexibleSECylinder_1*>(db.surfaces[surf2_ID - 1]);
	aA = &surf1->a;
	aB = &surf2->a;
	bA = &surf1->b;
	bB = &surf2->b;
	eA = &surf1->e;
	eB = &surf2->e;
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_A->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_A->getMatrix();
	dduiB = surf2->ddui_A->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_AAi->getMatrix();
	xBBi = surf2->x_BAi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	QBBi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
		QBBi[i] = new double[3];
	}
	//Salvando variáveis locais para montagem de superfícies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_AAi->MatrixToPtr(QABi, 3);
	surf2->Q_BAi->MatrixToPtr(QBBi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	double* ci = copy_convective;
	double* cp = convective;
	double* a4;
	double* a5;
	double* a6;
	double value = 0.0;
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

#pragma region AceGen
	double v01; double v010; double v011; double v012; double v013; double v014;
	double v015; double v016; double v017; double v018; double v019; double v02;
	double v020; double v021; double v022; double v023; double v024; double v025;
	double v026; double v027; double v03; double v04; double v05; double v06; double v07;
	double v08; double v09;
	int i1663, i2234, i4002, i4122, i4980, i4981, i4982, i4983, i4984, i4985, i4994, i4995
		, i4996, i4997, i4998, i4999, b38, b39, b1307, b1308, b1340, b1495, b1496, b1497, b1513
		, b1530, b1531, b1543, b1668, b1669, b1670, b1780, b1829, b2682, b2728, b2839, b2840, b2856
		, b3370, b3407;
	v[1] = cp[0];
	v[2] = cp[1];
	v[1564] = cos(v[2]);
	v[1457] = sin(v[2]);
	v[1458] = _copysign(1.e0, v[1457]);
	v[1459] = v[1458] * v[1564];
	v[506] = fabs(v[1457]);
	v[1460] = _copysign(1.e0, v[1564]);
	v[1461] = -(v[1457] * v[1460]);
	v[504] = fabs(v[1564]);
	v[3] = cp[2];
	v[4] = cp[3];
	v[1565] = cos(v[4]);
	v[1462] = sin(v[4]);
	v[1463] = _copysign(1.e0, v[1462]);
	v[1464] = v[1463] * v[1565];
	v[534] = fabs(v[1462]);
	v[1465] = _copysign(1.e0, v[1565]);
	v[1466] = -(v[1462] * v[1465]);
	v[532] = fabs(v[1565]);
	v[5] = ci[0];
	v[6] = ci[1];
	v[245] = sin(v[6]);
	v[243] = cos(v[6]);
	v[7] = ci[2];
	v[8] = ci[3];
	v[332] = sin(v[8]);
	v[330] = cos(v[8]);
	v[9] = gti[0];
	v[10] = gti[1];
	v[11] = gti[2];
	v[12] = (*epsn);
	v[13] = (*epsn0);
	v[14] = (*epst);
	v[15] = (*mus);
	v[16] = (*mud);
	v[17] = (*cn);
	v[18] = (*ct);
	v[4896] = v[18] / 2e0;
	v[19] = (*a4);
	v[4971] = v[18] * v[19];
	v[20] = (*a5);
	v[21] = (*a6);
	v[22] = invH[0][0];
	v[23] = invH[0][1];
	v[24] = invH[0][2];
	v[25] = invH[0][3];
	v[26] = invH[1][0];
	v[27] = invH[1][1];
	v[28] = invH[1][2];
	v[29] = invH[1][3];
	v[30] = invH[2][0];
	v[31] = invH[2][1];
	v[32] = invH[2][2];
	v[33] = invH[2][3];
	v[34] = invH[3][0];
	v[35] = invH[3][1];
	v[36] = invH[3][2];
	v[37] = invH[3][3];
	b38 = (*stick);
	b39 = (*previouscontact);
	v[40] = dA[0];
	v[41] = dA[1];
	v[42] = dA[2];
	v[43] = dA[3];
	v[986] = 2e0 * v[43];
	v[982] = v[43] / 2e0;
	v[183] = (v[43] * v[43]);
	v[44] = dA[4];
	v[985] = 2e0 * v[44];
	v[6662] = 0e0;
	v[6663] = 0e0;
	v[6664] = 0e0;
	v[6665] = -v[986];
	v[6666] = -v[985];
	v[6667] = 0e0;
	v[6668] = 0e0;
	v[6669] = 0e0;
	v[6670] = 0e0;
	v[6671] = 0e0;
	v[6672] = 0e0;
	v[6673] = 0e0;
	v[6674] = 0e0;
	v[6675] = 0e0;
	v[6676] = 0e0;
	v[6677] = 0e0;
	v[6678] = 0e0;
	v[6679] = 0e0;
	v[6680] = 0e0;
	v[6681] = 0e0;
	v[6682] = 0e0;
	v[6683] = 0e0;
	v[6684] = 0e0;
	v[6685] = 0e0;
	v[984] = v[44] / 2e0;
	v[6158] = 0e0;
	v[6159] = 0e0;
	v[6160] = 0e0;
	v[6161] = v[984];
	v[6162] = v[982];
	v[6163] = 0e0;
	v[6164] = 0e0;
	v[6165] = 0e0;
	v[6166] = 0e0;
	v[6167] = 0e0;
	v[6168] = 0e0;
	v[6169] = 0e0;
	v[6170] = 0e0;
	v[6171] = 0e0;
	v[6172] = 0e0;
	v[6173] = 0e0;
	v[6174] = 0e0;
	v[6175] = 0e0;
	v[6176] = 0e0;
	v[6177] = 0e0;
	v[6178] = 0e0;
	v[6179] = 0e0;
	v[6180] = 0e0;
	v[6181] = 0e0;
	v[181] = v[43] * v[984];
	v[176] = (v[44] * v[44]);
	v[932] = -v[176] - v[183];
	v[4799] = v[932] / 2e0;
	v[45] = dA[5];
	v[983] = v[45] / 2e0;
	v[6206] = 0e0;
	v[6207] = 0e0;
	v[6208] = 0e0;
	v[6209] = v[983];
	v[6210] = 0e0;
	v[6211] = v[982];
	v[6212] = 0e0;
	v[6213] = 0e0;
	v[6214] = 0e0;
	v[6215] = 0e0;
	v[6216] = 0e0;
	v[6217] = 0e0;
	v[6218] = 0e0;
	v[6219] = 0e0;
	v[6220] = 0e0;
	v[6221] = 0e0;
	v[6222] = 0e0;
	v[6223] = 0e0;
	v[6224] = 0e0;
	v[6225] = 0e0;
	v[6226] = 0e0;
	v[6227] = 0e0;
	v[6228] = 0e0;
	v[6229] = 0e0;
	v[6182] = 0e0;
	v[6183] = 0e0;
	v[6184] = 0e0;
	v[6185] = 0e0;
	v[6186] = v[983];
	v[6187] = v[984];
	v[6188] = 0e0;
	v[6189] = 0e0;
	v[6190] = 0e0;
	v[6191] = 0e0;
	v[6192] = 0e0;
	v[6193] = 0e0;
	v[6194] = 0e0;
	v[6195] = 0e0;
	v[6196] = 0e0;
	v[6197] = 0e0;
	v[6198] = 0e0;
	v[6199] = 0e0;
	v[6200] = 0e0;
	v[6201] = 0e0;
	v[6202] = 0e0;
	v[6203] = 0e0;
	v[6204] = 0e0;
	v[6205] = 0e0;
	v[981] = 2e0 * v[45];
	v[6566] = 0e0;
	v[6567] = 0e0;
	v[6568] = 0e0;
	v[6569] = -v[986];
	v[6570] = 0e0;
	v[6571] = -v[981];
	v[6572] = 0e0;
	v[6573] = 0e0;
	v[6574] = 0e0;
	v[6575] = 0e0;
	v[6576] = 0e0;
	v[6577] = 0e0;
	v[6578] = 0e0;
	v[6579] = 0e0;
	v[6580] = 0e0;
	v[6581] = 0e0;
	v[6582] = 0e0;
	v[6583] = 0e0;
	v[6584] = 0e0;
	v[6585] = 0e0;
	v[6586] = 0e0;
	v[6587] = 0e0;
	v[6588] = 0e0;
	v[6589] = 0e0;
	v[6470] = 0e0;
	v[6471] = 0e0;
	v[6472] = 0e0;
	v[6473] = 0e0;
	v[6474] = -v[985];
	v[6475] = -v[981];
	v[6476] = 0e0;
	v[6477] = 0e0;
	v[6478] = 0e0;
	v[6479] = 0e0;
	v[6480] = 0e0;
	v[6481] = 0e0;
	v[6482] = 0e0;
	v[6483] = 0e0;
	v[6484] = 0e0;
	v[6485] = 0e0;
	v[6486] = 0e0;
	v[6487] = 0e0;
	v[6488] = 0e0;
	v[6489] = 0e0;
	v[6490] = 0e0;
	v[6491] = 0e0;
	v[6492] = 0e0;
	v[6493] = 0e0;
	v[6446] = 0e0;
	v[6447] = 0e0;
	v[6448] = 0e0;
	v[6449] = v[986];
	v[6450] = v[985];
	v[6451] = v[981];
	v[6452] = 0e0;
	v[6453] = 0e0;
	v[6454] = 0e0;
	v[6455] = 0e0;
	v[6456] = 0e0;
	v[6457] = 0e0;
	v[6458] = 0e0;
	v[6459] = 0e0;
	v[6460] = 0e0;
	v[6461] = 0e0;
	v[6462] = 0e0;
	v[6463] = 0e0;
	v[6464] = 0e0;
	v[6465] = 0e0;
	v[6466] = 0e0;
	v[6467] = 0e0;
	v[6468] = 0e0;
	v[6469] = 0e0;
	v[950] = v[181] - v[45];
	v[944] = v[181] + v[45];
	v[188] = v[44] * v[983];
	v[937] = v[188] - v[43];
	v[935] = v[188] + v[43];
	v[186] = v[43] * v[983];
	v[945] = v[186] + v[44];
	v[936] = v[186] - v[44];
	v[177] = (v[45] * v[45]);
	v[961] = 4e0 + v[176] + v[177] + v[183];
	v[5311] = 8e0 / Power(v[961], 3);
	v[962] = -4e0 / (v[961] * v[961]);
	v[955] = -v[176] - v[177];
	v[4801] = v[955] / 2e0;
	v[942] = -v[177] - v[183];
	v[4800] = v[942] / 2e0;
	v[46] = dA[6];
	v[47] = dA[7];
	v[48] = dA[8];
	v[49] = dA[9];
	v[980] = 2e0 * v[49];
	v[976] = v[49] / 2e0;
	v[202] = (v[49] * v[49]);
	v[50] = dA[10];
	v[979] = 2e0 * v[50];
	v[6974] = 0e0;
	v[6975] = 0e0;
	v[6976] = 0e0;
	v[6977] = 0e0;
	v[6978] = 0e0;
	v[6979] = 0e0;
	v[6980] = 0e0;
	v[6981] = 0e0;
	v[6982] = 0e0;
	v[6983] = -v[980];
	v[6984] = -v[979];
	v[6985] = 0e0;
	v[6986] = 0e0;
	v[6987] = 0e0;
	v[6988] = 0e0;
	v[6989] = 0e0;
	v[6990] = 0e0;
	v[6991] = 0e0;
	v[6992] = 0e0;
	v[6993] = 0e0;
	v[6994] = 0e0;
	v[6995] = 0e0;
	v[6996] = 0e0;
	v[6997] = 0e0;
	v[978] = v[50] / 2e0;
	v[6230] = 0e0;
	v[6231] = 0e0;
	v[6232] = 0e0;
	v[6233] = 0e0;
	v[6234] = 0e0;
	v[6235] = 0e0;
	v[6236] = 0e0;
	v[6237] = 0e0;
	v[6238] = 0e0;
	v[6239] = v[978];
	v[6240] = v[976];
	v[6241] = 0e0;
	v[6242] = 0e0;
	v[6243] = 0e0;
	v[6244] = 0e0;
	v[6245] = 0e0;
	v[6246] = 0e0;
	v[6247] = 0e0;
	v[6248] = 0e0;
	v[6249] = 0e0;
	v[6250] = 0e0;
	v[6251] = 0e0;
	v[6252] = 0e0;
	v[6253] = 0e0;
	v[200] = v[49] * v[978];
	v[195] = (v[50] * v[50]);
	v[901] = -v[195] - v[202];
	v[4796] = v[901] / 2e0;
	v[51] = dA[11];
	v[977] = v[51] / 2e0;
	v[6278] = 0e0;
	v[6279] = 0e0;
	v[6280] = 0e0;
	v[6281] = 0e0;
	v[6282] = 0e0;
	v[6283] = 0e0;
	v[6284] = 0e0;
	v[6285] = 0e0;
	v[6286] = 0e0;
	v[6287] = v[977];
	v[6288] = 0e0;
	v[6289] = v[976];
	v[6290] = 0e0;
	v[6291] = 0e0;
	v[6292] = 0e0;
	v[6293] = 0e0;
	v[6294] = 0e0;
	v[6295] = 0e0;
	v[6296] = 0e0;
	v[6297] = 0e0;
	v[6298] = 0e0;
	v[6299] = 0e0;
	v[6300] = 0e0;
	v[6301] = 0e0;
	v[6254] = 0e0;
	v[6255] = 0e0;
	v[6256] = 0e0;
	v[6257] = 0e0;
	v[6258] = 0e0;
	v[6259] = 0e0;
	v[6260] = 0e0;
	v[6261] = 0e0;
	v[6262] = 0e0;
	v[6263] = 0e0;
	v[6264] = v[977];
	v[6265] = v[978];
	v[6266] = 0e0;
	v[6267] = 0e0;
	v[6268] = 0e0;
	v[6269] = 0e0;
	v[6270] = 0e0;
	v[6271] = 0e0;
	v[6272] = 0e0;
	v[6273] = 0e0;
	v[6274] = 0e0;
	v[6275] = 0e0;
	v[6276] = 0e0;
	v[6277] = 0e0;
	v[975] = 2e0 * v[51];
	v[6878] = 0e0;
	v[6879] = 0e0;
	v[6880] = 0e0;
	v[6881] = 0e0;
	v[6882] = 0e0;
	v[6883] = 0e0;
	v[6884] = 0e0;
	v[6885] = 0e0;
	v[6886] = 0e0;
	v[6887] = -v[980];
	v[6888] = 0e0;
	v[6889] = -v[975];
	v[6890] = 0e0;
	v[6891] = 0e0;
	v[6892] = 0e0;
	v[6893] = 0e0;
	v[6894] = 0e0;
	v[6895] = 0e0;
	v[6896] = 0e0;
	v[6897] = 0e0;
	v[6898] = 0e0;
	v[6899] = 0e0;
	v[6900] = 0e0;
	v[6901] = 0e0;
	v[6782] = 0e0;
	v[6783] = 0e0;
	v[6784] = 0e0;
	v[6785] = 0e0;
	v[6786] = 0e0;
	v[6787] = 0e0;
	v[6788] = 0e0;
	v[6789] = 0e0;
	v[6790] = 0e0;
	v[6791] = 0e0;
	v[6792] = -v[979];
	v[6793] = -v[975];
	v[6794] = 0e0;
	v[6795] = 0e0;
	v[6796] = 0e0;
	v[6797] = 0e0;
	v[6798] = 0e0;
	v[6799] = 0e0;
	v[6800] = 0e0;
	v[6801] = 0e0;
	v[6802] = 0e0;
	v[6803] = 0e0;
	v[6804] = 0e0;
	v[6805] = 0e0;
	v[6758] = 0e0;
	v[6759] = 0e0;
	v[6760] = 0e0;
	v[6761] = 0e0;
	v[6762] = 0e0;
	v[6763] = 0e0;
	v[6764] = 0e0;
	v[6765] = 0e0;
	v[6766] = 0e0;
	v[6767] = v[980];
	v[6768] = v[979];
	v[6769] = v[975];
	v[6770] = 0e0;
	v[6771] = 0e0;
	v[6772] = 0e0;
	v[6773] = 0e0;
	v[6774] = 0e0;
	v[6775] = 0e0;
	v[6776] = 0e0;
	v[6777] = 0e0;
	v[6778] = 0e0;
	v[6779] = 0e0;
	v[6780] = 0e0;
	v[6781] = 0e0;
	v[919] = v[200] - v[51];
	v[913] = v[200] + v[51];
	v[207] = v[50] * v[977];
	v[906] = v[207] - v[49];
	v[904] = v[207] + v[49];
	v[205] = v[49] * v[977];
	v[914] = v[205] + v[50];
	v[905] = v[205] - v[50];
	v[196] = (v[51] * v[51]);
	v[930] = 4e0 + v[195] + v[196] + v[202];
	v[5307] = 8e0 / Power(v[930], 3);
	v[931] = -4e0 / (v[930] * v[930]);
	v[924] = -v[195] - v[196];
	v[4798] = v[924] / 2e0;
	v[911] = -v[196] - v[202];
	v[4797] = v[911] / 2e0;
	v[367] = duiA[3] * v[20] + dduiA[3] * v[21] + v[19] * v[43];
	v[373] = duiA[4] * v[20] + dduiA[4] * v[21] + v[19] * v[44];
	v[375] = duiA[5] * v[20] + dduiA[5] * v[21] + v[19] * v[45];
	v[393] = duiA[9] * v[20] + dduiA[9] * v[21] + v[19] * v[49];
	v[399] = duiA[10] * v[20] + dduiA[10] * v[21] + v[19] * v[50];
	v[401] = duiA[11] * v[20] + dduiA[11] * v[21] + v[19] * v[51];
	v[76] = (*aA);
	v[5105] = v[1564] * v[76];
	v[77] = (*bA);
	v[4735] = v[1457] * v[77];
	v[79] = xAAi[0];
	v[80] = xAAi[1];
	v[81] = xAAi[2];
	v[82] = xBAi[0];
	v[83] = xBAi[1];
	v[84] = xBAi[2];
	v[85] = QAAi[0][0];
	v[86] = QAAi[0][1];
	v[88] = QAAi[1][0];
	v[89] = QAAi[1][1];
	v[91] = QAAi[2][0];
	v[92] = QAAi[2][1];
	v[94] = QBAi[0][0];
	v[95] = QBAi[0][1];
	v[97] = QBAi[1][0];
	v[98] = QBAi[1][1];
	v[100] = QBAi[2][0];
	v[101] = QBAi[2][1];
	v[103] = dB[0];
	v[104] = dB[1];
	v[105] = dB[2];
	v[106] = dB[3];
	v[974] = 2e0 * v[106];
	v[970] = v[106] / 2e0;
	v[270] = (v[106] * v[106]);
	v[107] = dB[4];
	v[973] = 2e0 * v[107];
	v[7430] = 0e0;
	v[7431] = 0e0;
	v[7432] = 0e0;
	v[7433] = 0e0;
	v[7434] = 0e0;
	v[7435] = 0e0;
	v[7436] = 0e0;
	v[7437] = 0e0;
	v[7438] = 0e0;
	v[7439] = 0e0;
	v[7440] = 0e0;
	v[7441] = 0e0;
	v[7442] = 0e0;
	v[7443] = 0e0;
	v[7444] = 0e0;
	v[7445] = -v[974];
	v[7446] = -v[973];
	v[7447] = 0e0;
	v[7448] = 0e0;
	v[7449] = 0e0;
	v[7450] = 0e0;
	v[7451] = 0e0;
	v[7452] = 0e0;
	v[7453] = 0e0;
	v[972] = v[107] / 2e0;
	v[6302] = 0e0;
	v[6303] = 0e0;
	v[6304] = 0e0;
	v[6305] = 0e0;
	v[6306] = 0e0;
	v[6307] = 0e0;
	v[6308] = 0e0;
	v[6309] = 0e0;
	v[6310] = 0e0;
	v[6311] = 0e0;
	v[6312] = 0e0;
	v[6313] = 0e0;
	v[6314] = 0e0;
	v[6315] = 0e0;
	v[6316] = 0e0;
	v[6317] = v[972];
	v[6318] = v[970];
	v[6319] = 0e0;
	v[6320] = 0e0;
	v[6321] = 0e0;
	v[6322] = 0e0;
	v[6323] = 0e0;
	v[6324] = 0e0;
	v[6325] = 0e0;
	v[268] = v[106] * v[972];
	v[263] = (v[107] * v[107]);
	v[750] = -v[263] - v[270];
	v[4787] = v[750] / 2e0;
	v[108] = dB[5];
	v[971] = v[108] / 2e0;
	v[6350] = 0e0;
	v[6351] = 0e0;
	v[6352] = 0e0;
	v[6353] = 0e0;
	v[6354] = 0e0;
	v[6355] = 0e0;
	v[6356] = 0e0;
	v[6357] = 0e0;
	v[6358] = 0e0;
	v[6359] = 0e0;
	v[6360] = 0e0;
	v[6361] = 0e0;
	v[6362] = 0e0;
	v[6363] = 0e0;
	v[6364] = 0e0;
	v[6365] = v[971];
	v[6366] = 0e0;
	v[6367] = v[970];
	v[6368] = 0e0;
	v[6369] = 0e0;
	v[6370] = 0e0;
	v[6371] = 0e0;
	v[6372] = 0e0;
	v[6373] = 0e0;
	v[6326] = 0e0;
	v[6327] = 0e0;
	v[6328] = 0e0;
	v[6329] = 0e0;
	v[6330] = 0e0;
	v[6331] = 0e0;
	v[6332] = 0e0;
	v[6333] = 0e0;
	v[6334] = 0e0;
	v[6335] = 0e0;
	v[6336] = 0e0;
	v[6337] = 0e0;
	v[6338] = 0e0;
	v[6339] = 0e0;
	v[6340] = 0e0;
	v[6341] = 0e0;
	v[6342] = v[971];
	v[6343] = v[972];
	v[6344] = 0e0;
	v[6345] = 0e0;
	v[6346] = 0e0;
	v[6347] = 0e0;
	v[6348] = 0e0;
	v[6349] = 0e0;
	v[969] = 2e0 * v[108];
	v[7334] = 0e0;
	v[7335] = 0e0;
	v[7336] = 0e0;
	v[7337] = 0e0;
	v[7338] = 0e0;
	v[7339] = 0e0;
	v[7340] = 0e0;
	v[7341] = 0e0;
	v[7342] = 0e0;
	v[7343] = 0e0;
	v[7344] = 0e0;
	v[7345] = 0e0;
	v[7346] = 0e0;
	v[7347] = 0e0;
	v[7348] = 0e0;
	v[7349] = -v[974];
	v[7350] = 0e0;
	v[7351] = -v[969];
	v[7352] = 0e0;
	v[7353] = 0e0;
	v[7354] = 0e0;
	v[7355] = 0e0;
	v[7356] = 0e0;
	v[7357] = 0e0;
	v[7238] = 0e0;
	v[7239] = 0e0;
	v[7240] = 0e0;
	v[7241] = 0e0;
	v[7242] = 0e0;
	v[7243] = 0e0;
	v[7244] = 0e0;
	v[7245] = 0e0;
	v[7246] = 0e0;
	v[7247] = 0e0;
	v[7248] = 0e0;
	v[7249] = 0e0;
	v[7250] = 0e0;
	v[7251] = 0e0;
	v[7252] = 0e0;
	v[7253] = 0e0;
	v[7254] = -v[973];
	v[7255] = -v[969];
	v[7256] = 0e0;
	v[7257] = 0e0;
	v[7258] = 0e0;
	v[7259] = 0e0;
	v[7260] = 0e0;
	v[7261] = 0e0;
	v[7214] = 0e0;
	v[7215] = 0e0;
	v[7216] = 0e0;
	v[7217] = 0e0;
	v[7218] = 0e0;
	v[7219] = 0e0;
	v[7220] = 0e0;
	v[7221] = 0e0;
	v[7222] = 0e0;
	v[7223] = 0e0;
	v[7224] = 0e0;
	v[7225] = 0e0;
	v[7226] = 0e0;
	v[7227] = 0e0;
	v[7228] = 0e0;
	v[7229] = v[974];
	v[7230] = v[973];
	v[7231] = v[969];
	v[7232] = 0e0;
	v[7233] = 0e0;
	v[7234] = 0e0;
	v[7235] = 0e0;
	v[7236] = 0e0;
	v[7237] = 0e0;
	v[768] = -v[108] + v[268];
	v[762] = v[108] + v[268];
	v[275] = v[107] * v[971];
	v[755] = -v[106] + v[275];
	v[753] = v[106] + v[275];
	v[273] = v[106] * v[971];
	v[763] = v[107] + v[273];
	v[754] = -v[107] + v[273];
	v[264] = (v[108] * v[108]);
	v[779] = 4e0 + v[263] + v[264] + v[270];
	v[5303] = 8e0 / Power(v[779], 3);
	v[780] = -4e0 / (v[779] * v[779]);
	v[773] = -v[263] - v[264];
	v[4789] = v[773] / 2e0;
	v[760] = -v[264] - v[270];
	v[4788] = v[760] / 2e0;
	v[109] = dB[6];
	v[110] = dB[7];
	v[111] = dB[8];
	v[112] = dB[9];
	v[968] = 2e0 * v[112];
	v[964] = v[112] / 2e0;
	v[289] = (v[112] * v[112]);
	v[113] = dB[10];
	v[967] = 2e0 * v[113];
	v[7742] = 0e0;
	v[7743] = 0e0;
	v[7744] = 0e0;
	v[7745] = 0e0;
	v[7746] = 0e0;
	v[7747] = 0e0;
	v[7748] = 0e0;
	v[7749] = 0e0;
	v[7750] = 0e0;
	v[7751] = 0e0;
	v[7752] = 0e0;
	v[7753] = 0e0;
	v[7754] = 0e0;
	v[7755] = 0e0;
	v[7756] = 0e0;
	v[7757] = 0e0;
	v[7758] = 0e0;
	v[7759] = 0e0;
	v[7760] = 0e0;
	v[7761] = 0e0;
	v[7762] = 0e0;
	v[7763] = -v[968];
	v[7764] = -v[967];
	v[7765] = 0e0;
	v[966] = v[113] / 2e0;
	v[6374] = 0e0;
	v[6375] = 0e0;
	v[6376] = 0e0;
	v[6377] = 0e0;
	v[6378] = 0e0;
	v[6379] = 0e0;
	v[6380] = 0e0;
	v[6381] = 0e0;
	v[6382] = 0e0;
	v[6383] = 0e0;
	v[6384] = 0e0;
	v[6385] = 0e0;
	v[6386] = 0e0;
	v[6387] = 0e0;
	v[6388] = 0e0;
	v[6389] = 0e0;
	v[6390] = 0e0;
	v[6391] = 0e0;
	v[6392] = 0e0;
	v[6393] = 0e0;
	v[6394] = 0e0;
	v[6395] = v[966];
	v[6396] = v[964];
	v[6397] = 0e0;
	v[287] = v[112] * v[966];
	v[282] = (v[113] * v[113]);
	v[719] = -v[282] - v[289];
	v[4784] = v[719] / 2e0;
	v[114] = dB[11];
	v[965] = v[114] / 2e0;
	v[6422] = 0e0;
	v[6423] = 0e0;
	v[6424] = 0e0;
	v[6425] = 0e0;
	v[6426] = 0e0;
	v[6427] = 0e0;
	v[6428] = 0e0;
	v[6429] = 0e0;
	v[6430] = 0e0;
	v[6431] = 0e0;
	v[6432] = 0e0;
	v[6433] = 0e0;
	v[6434] = 0e0;
	v[6435] = 0e0;
	v[6436] = 0e0;
	v[6437] = 0e0;
	v[6438] = 0e0;
	v[6439] = 0e0;
	v[6440] = 0e0;
	v[6441] = 0e0;
	v[6442] = 0e0;
	v[6443] = v[965];
	v[6444] = 0e0;
	v[6445] = v[964];
	v[6398] = 0e0;
	v[6399] = 0e0;
	v[6400] = 0e0;
	v[6401] = 0e0;
	v[6402] = 0e0;
	v[6403] = 0e0;
	v[6404] = 0e0;
	v[6405] = 0e0;
	v[6406] = 0e0;
	v[6407] = 0e0;
	v[6408] = 0e0;
	v[6409] = 0e0;
	v[6410] = 0e0;
	v[6411] = 0e0;
	v[6412] = 0e0;
	v[6413] = 0e0;
	v[6414] = 0e0;
	v[6415] = 0e0;
	v[6416] = 0e0;
	v[6417] = 0e0;
	v[6418] = 0e0;
	v[6419] = 0e0;
	v[6420] = v[965];
	v[6421] = v[966];
	v[963] = 2e0 * v[114];
	v[7646] = 0e0;
	v[7647] = 0e0;
	v[7648] = 0e0;
	v[7649] = 0e0;
	v[7650] = 0e0;
	v[7651] = 0e0;
	v[7652] = 0e0;
	v[7653] = 0e0;
	v[7654] = 0e0;
	v[7655] = 0e0;
	v[7656] = 0e0;
	v[7657] = 0e0;
	v[7658] = 0e0;
	v[7659] = 0e0;
	v[7660] = 0e0;
	v[7661] = 0e0;
	v[7662] = 0e0;
	v[7663] = 0e0;
	v[7664] = 0e0;
	v[7665] = 0e0;
	v[7666] = 0e0;
	v[7667] = -v[968];
	v[7668] = 0e0;
	v[7669] = -v[963];
	v[7550] = 0e0;
	v[7551] = 0e0;
	v[7552] = 0e0;
	v[7553] = 0e0;
	v[7554] = 0e0;
	v[7555] = 0e0;
	v[7556] = 0e0;
	v[7557] = 0e0;
	v[7558] = 0e0;
	v[7559] = 0e0;
	v[7560] = 0e0;
	v[7561] = 0e0;
	v[7562] = 0e0;
	v[7563] = 0e0;
	v[7564] = 0e0;
	v[7565] = 0e0;
	v[7566] = 0e0;
	v[7567] = 0e0;
	v[7568] = 0e0;
	v[7569] = 0e0;
	v[7570] = 0e0;
	v[7571] = 0e0;
	v[7572] = -v[967];
	v[7573] = -v[963];
	v[7526] = 0e0;
	v[7527] = 0e0;
	v[7528] = 0e0;
	v[7529] = 0e0;
	v[7530] = 0e0;
	v[7531] = 0e0;
	v[7532] = 0e0;
	v[7533] = 0e0;
	v[7534] = 0e0;
	v[7535] = 0e0;
	v[7536] = 0e0;
	v[7537] = 0e0;
	v[7538] = 0e0;
	v[7539] = 0e0;
	v[7540] = 0e0;
	v[7541] = 0e0;
	v[7542] = 0e0;
	v[7543] = 0e0;
	v[7544] = 0e0;
	v[7545] = 0e0;
	v[7546] = 0e0;
	v[7547] = v[968];
	v[7548] = v[967];
	v[7549] = v[963];
	v[737] = -v[114] + v[287];
	v[731] = v[114] + v[287];
	v[294] = v[113] * v[965];
	v[724] = -v[112] + v[294];
	v[722] = v[112] + v[294];
	v[292] = v[112] * v[965];
	v[732] = v[113] + v[292];
	v[723] = -v[113] + v[292];
	v[283] = (v[114] * v[114]);
	v[748] = 4e0 + v[282] + v[283] + v[289];
	v[5299] = 8e0 / Power(v[748], 3);
	v[749] = -4e0 / (v[748] * v[748]);
	v[742] = -v[282] - v[283];
	v[4786] = v[742] / 2e0;
	v[729] = -v[283] - v[289];
	v[4785] = v[729] / 2e0;
	v[419] = v[106] * v[19] + duiB[3] * v[20] + dduiB[3] * v[21];
	v[425] = v[107] * v[19] + duiB[4] * v[20] + dduiB[4] * v[21];
	v[427] = v[108] * v[19] + duiB[5] * v[20] + dduiB[5] * v[21];
	v[445] = v[112] * v[19] + duiB[9] * v[20] + dduiB[9] * v[21];
	v[451] = v[113] * v[19] + duiB[10] * v[20] + dduiB[10] * v[21];
	v[453] = v[114] * v[19] + duiB[11] * v[20] + dduiB[11] * v[21];
	v[139] = (*aB);
	v[5087] = v[139] * v[1565];
	v[140] = (*bB);
	v[4748] = v[140] * v[1462];
	v[142] = xABi[0];
	v[143] = xABi[1];
	v[144] = xABi[2];
	v[145] = xBBi[0];
	v[146] = xBBi[1];
	v[147] = xBBi[2];
	v[148] = QABi[0][0];
	v[149] = QABi[0][1];
	v[151] = QABi[1][0];
	v[152] = QABi[1][1];
	v[154] = QABi[2][0];
	v[155] = QABi[2][1];
	v[157] = QBBi[0][0];
	v[158] = QBBi[0][1];
	v[160] = QBBi[1][0];
	v[161] = QBBi[1][1];
	v[163] = QBBi[2][0];
	v[164] = QBBi[2][1];
	v[175] = 4e0 / v[961];
	v[5152] = 2e0 * v[175];
	v[4724] = v[175] / 2e0;
	v[364] = (v[175] * v[175]);
	v[5458] = 1e0 / Power(v[364], 3);
	v[4726] = v[375] / v[364];
	v[4725] = v[373] / v[364];
	v[4723] = v[367] / v[364];
	v[178] = 1e0 + v[4724] * v[955];
	v[1986] = v[178] * v[4723];
	v[179] = v[175] * v[950];
	v[1988] = v[179] * v[4725];
	v[2770] = v[1986] + v[1988];
	v[180] = v[175] * v[945];
	v[5256] = v[180] / v[364];
	v[1989] = v[180] * v[4726];
	v[2775] = v[1986] + v[1989];
	v[2764] = v[1988] + v[2775];
	v[182] = v[175] * v[944];
	v[1993] = v[182] * v[4723];
	v[184] = 1e0 + v[4724] * v[942];
	v[1981] = v[184] * v[4725];
	v[2766] = v[1981] + v[1993];
	v[185] = v[175] * v[937];
	v[5257] = v[185] / v[364];
	v[1982] = v[185] * v[4726];
	v[2776] = v[1981] + v[1982];
	v[2769] = v[1993] + v[2776];
	v[187] = v[175] * v[936];
	v[1997] = v[187] * v[4723];
	v[189] = v[175] * v[935];
	v[5258] = v[189] / v[364];
	v[1977] = v[189] * v[4725];
	v[190] = 1e0 + v[4724] * v[932];
	v[1978] = v[190] * v[4726];
	v[2774] = v[1977] + v[1978] + v[1997];
	v[2771] = -v[1997] + v[2774];
	v[2765] = -v[1977] + v[2774];
	v[191] = -(v[45] * v[4724]);
	v[4973] = 2e0 * v[191];
	v[370] = -(v[175] * v[191]);
	v[192] = v[44] * v[4724];
	v[4974] = 2e0 * v[192];
	v[382] = -(v[191] * v[192]);
	v[379] = -(v[175] * v[192]);
	v[193] = -(v[43] * v[4724]);
	v[4975] = 2e0 * v[193];
	v[383] = -(v[175] * v[193]);
	v[4753] = v[382] - v[383];
	v[380] = v[191] * v[193];
	v[4751] = -v[379] + v[380];
	v[371] = -(v[192] * v[193]);
	v[4749] = -v[370] + v[371];
	v[194] = 4e0 / v[930];
	v[5139] = 2e0 * v[194];
	v[4728] = v[194] / 2e0;
	v[390] = (v[194] * v[194]);
	v[5455] = 1e0 / Power(v[390], 3);
	v[4730] = v[401] / v[390];
	v[4729] = v[399] / v[390];
	v[4727] = v[393] / v[390];
	v[197] = 1e0 + v[4728] * v[924];
	v[1959] = v[197] * v[4727];
	v[198] = v[194] * v[919];
	v[1961] = v[198] * v[4729];
	v[2785] = v[1959] + v[1961];
	v[199] = v[194] * v[914];
	v[5253] = v[199] / v[390];
	v[1962] = v[199] * v[4730];
	v[2790] = v[1959] + v[1962];
	v[2779] = v[1961] + v[2790];
	v[201] = v[194] * v[913];
	v[1966] = v[201] * v[4727];
	v[203] = 1e0 + v[4728] * v[911];
	v[1954] = v[203] * v[4729];
	v[2781] = v[1954] + v[1966];
	v[204] = v[194] * v[906];
	v[5254] = v[204] / v[390];
	v[1955] = v[204] * v[4730];
	v[2791] = v[1954] + v[1955];
	v[2784] = v[1966] + v[2791];
	v[206] = v[194] * v[905];
	v[1970] = v[206] * v[4727];
	v[208] = v[194] * v[904];
	v[5255] = v[208] / v[390];
	v[1950] = v[208] * v[4729];
	v[209] = 1e0 + v[4728] * v[901];
	v[1951] = v[209] * v[4730];
	v[2789] = v[1950] + v[1951] + v[1970];
	v[2786] = -v[1970] + v[2789];
	v[2780] = -v[1950] + v[2789];
	v[210] = -(v[4728] * v[51]);
	v[4977] = 2e0 * v[210];
	v[396] = -(v[194] * v[210]);
	v[211] = v[4728] * v[50];
	v[4978] = 2e0 * v[211];
	v[408] = -(v[210] * v[211]);
	v[405] = -(v[194] * v[211]);
	v[212] = -(v[4728] * v[49]);
	v[4979] = 2e0 * v[212];
	v[409] = -(v[194] * v[212]);
	v[4758] = v[408] - v[409];
	v[406] = v[210] * v[212];
	v[4756] = -v[405] + v[406];
	v[397] = -(v[211] * v[212]);
	v[4754] = -v[396] + v[397];
	v[213] = v[178] * v[85] + v[179] * v[88] + v[180] * v[91];
	v[214] = v[178] * v[86] + v[179] * v[89] + v[180] * v[92];
	v[216] = v[182] * v[85] + v[184] * v[88] + v[185] * v[91];
	v[217] = v[182] * v[86] + v[184] * v[89] + v[185] * v[92];
	v[219] = v[187] * v[85] + v[189] * v[88] + v[190] * v[91];
	v[220] = v[187] * v[86] + v[189] * v[89] + v[190] * v[92];
	v[222] = v[100] * v[199] + v[197] * v[94] + v[198] * v[97];
	v[223] = v[101] * v[199] + v[197] * v[95] + v[198] * v[98];
	v[225] = v[100] * v[204] + v[201] * v[94] + v[203] * v[97];
	v[226] = v[101] * v[204] + v[201] * v[95] + v[203] * v[98];
	v[228] = v[100] * v[209] + v[206] * v[94] + v[208] * v[97];
	v[229] = v[101] * v[209] + v[206] * v[95] + v[208] * v[98];
	v[231] = v[40] + v[79];
	v[232] = v[41] + v[80];
	v[233] = v[42] + v[81];
	v[234] = v[46] + v[82];
	v[235] = v[47] + v[83];
	v[236] = v[48] + v[84];
	v[237] = (1e0 - v[5]) / 2e0;
	v[238] = (1e0 + v[5]) / 2e0;
	v[239] = (1e0 - v[1]) / 2e0;
	v[240] = (1e0 + v[1]) / 2e0;
	v[4048] = v[219] * v[239] + v[228] * v[240];
	v[4047] = v[216] * v[239] + v[225] * v[240];
	v[4046] = v[213] * v[239] + v[222] * v[240];
	v[4044] = v[220] * v[239] + v[229] * v[240];
	v[4043] = v[217] * v[239] + v[226] * v[240];
	v[4042] = v[214] * v[239] + v[223] * v[240];
	v[241] = 2e0 / (*eA);
	v[4731] = -1e0 + v[241];
	v[3973] = Power(v[504], v[4731]);
	v[5227] = -(v[1460] * v[3973]);
	v[3972] = -2e0 + v[241];
	v[3971] = Power(v[506], v[4731]);
	v[5230] = v[1458] * v[3971];
	v[1467] = v[241] * (v[1459] * v[3971] + v[1461] * v[3973]);
	v[4705] = v[1461] * v[4731] * Power(v[504], v[3972]);
	v[4703] = v[1459] * v[4731] * Power(v[506], v[3972]);
	v[503] = Power(v[504], v[241]) + Power(v[506], v[241]);
	v[3758] = -1e0 - 1e0 / v[241];
	v[5286] = -((v[1467] * v[3758] * Power(v[503], -1e0 + v[3758])) / v[241]);
	v[3757] = Power(v[503], v[3758]);
	v[5106] = -(v[3757] / v[241]);
	v[507] = v[1467] * v[5106];
	v[242] = 1e0 / Power(Power(fabs(v[243]), v[241]) + Power(fabs(v[245]), v[241]), 1 / v[241]);
	v[244] = v[242] * v[243] * v[76];
	v[246] = v[242] * v[245] * v[77];
	v[248] = 1e0 / Power(v[503], 1 / v[241]);
	v[5228] = v[248] * v[77];
	v[4734] = v[248] * v[76];
	v[4240] = v[1564] * v[5228];
	v[4238] = -(v[1457] * v[4734]);
	v[509] = v[4240] + v[4735] * v[507];
	v[508] = v[4238] + v[507] * v[5105];
	v[515] = v[213] * v[508] + v[214] * v[509];
	v[514] = v[222] * v[508] + v[223] * v[509];
	v[516] = v[240] * v[514] + v[239] * v[515];
	v[513] = v[216] * v[508] + v[217] * v[509];
	v[512] = v[225] * v[508] + v[226] * v[509];
	v[517] = v[240] * v[512] + v[239] * v[513];
	v[511] = v[219] * v[508] + v[220] * v[509];
	v[510] = v[228] * v[508] + v[229] * v[509];
	v[518] = v[240] * v[510] + v[239] * v[511];
	v[250] = v[1564] * v[4734];
	v[252] = v[248] * v[4735];
	v[3765] = v[100] * v[250] + v[101] * v[252];
	v[3764] = v[250] * v[97] + v[252] * v[98];
	v[3763] = v[250] * v[94] + v[252] * v[95];
	v[3762] = v[250] * v[91] + v[252] * v[92];
	v[3761] = v[250] * v[88] + v[252] * v[89];
	v[3760] = v[250] * v[85] + v[252] * v[86];
	v[499] = v[236] + v[228] * v[250] + v[229] * v[252];
	v[498] = v[233] + v[219] * v[250] + v[220] * v[252];
	v[500] = (-v[498] + v[499]) / 2e0;
	v[496] = v[235] + v[225] * v[250] + v[226] * v[252];
	v[495] = v[232] + v[216] * v[250] + v[217] * v[252];
	v[497] = (-v[495] + v[496]) / 2e0;
	v[493] = v[234] + v[222] * v[250] + v[223] * v[252];
	v[492] = v[231] + v[213] * v[250] + v[214] * v[252];
	v[494] = (-v[492] + v[493]) / 2e0;
	v[262] = 4e0 / v[779];
	v[5126] = 2e0 * v[262];
	v[4737] = v[262] / 2e0;
	v[416] = (v[262] * v[262]);
	v[5452] = 1e0 / Power(v[416], 3);
	v[4739] = v[427] / v[416];
	v[4738] = v[425] / v[416];
	v[4736] = v[419] / v[416];
	v[265] = 1e0 + v[4737] * v[773];
	v[1932] = v[265] * v[4736];
	v[266] = v[262] * v[768];
	v[1934] = v[266] * v[4738];
	v[2800] = v[1932] + v[1934];
	v[267] = v[262] * v[763];
	v[5250] = v[267] / v[416];
	v[1935] = v[267] * v[4739];
	v[2805] = v[1932] + v[1935];
	v[2794] = v[1934] + v[2805];
	v[269] = v[262] * v[762];
	v[1939] = v[269] * v[4736];
	v[271] = 1e0 + v[4737] * v[760];
	v[1927] = v[271] * v[4738];
	v[2796] = v[1927] + v[1939];
	v[272] = v[262] * v[755];
	v[5251] = v[272] / v[416];
	v[1928] = v[272] * v[4739];
	v[2806] = v[1927] + v[1928];
	v[2799] = v[1939] + v[2806];
	v[274] = v[262] * v[754];
	v[1943] = v[274] * v[4736];
	v[276] = v[262] * v[753];
	v[5252] = v[276] / v[416];
	v[1923] = v[276] * v[4738];
	v[277] = 1e0 + v[4737] * v[750];
	v[1924] = v[277] * v[4739];
	v[2804] = v[1923] + v[1924] + v[1943];
	v[2801] = -v[1943] + v[2804];
	v[2795] = -v[1923] + v[2804];
	v[278] = -(v[108] * v[4737]);
	v[4987] = 2e0 * v[278];
	v[422] = -(v[262] * v[278]);
	v[279] = v[107] * v[4737];
	v[4988] = 2e0 * v[279];
	v[434] = -(v[278] * v[279]);
	v[431] = -(v[262] * v[279]);
	v[280] = -(v[106] * v[4737]);
	v[4989] = 2e0 * v[280];
	v[435] = -(v[262] * v[280]);
	v[4763] = v[434] - v[435];
	v[432] = v[278] * v[280];
	v[4761] = -v[431] + v[432];
	v[423] = -(v[279] * v[280]);
	v[4759] = -v[422] + v[423];
	v[281] = 4e0 / v[748];
	v[5113] = 2e0 * v[281];
	v[4741] = v[281] / 2e0;
	v[442] = (v[281] * v[281]);
	v[5449] = 1e0 / Power(v[442], 3);
	v[4743] = v[453] / v[442];
	v[4742] = v[451] / v[442];
	v[4740] = v[445] / v[442];
	v[284] = 1e0 + v[4741] * v[742];
	v[1905] = v[284] * v[4740];
	v[285] = v[281] * v[737];
	v[1907] = v[285] * v[4742];
	v[2815] = v[1905] + v[1907];
	v[286] = v[281] * v[732];
	v[5247] = v[286] / v[442];
	v[1908] = v[286] * v[4743];
	v[2820] = v[1905] + v[1908];
	v[2809] = v[1907] + v[2820];
	v[288] = v[281] * v[731];
	v[1912] = v[288] * v[4740];
	v[290] = 1e0 + v[4741] * v[729];
	v[1900] = v[290] * v[4742];
	v[2811] = v[1900] + v[1912];
	v[291] = v[281] * v[724];
	v[5248] = v[291] / v[442];
	v[1901] = v[291] * v[4743];
	v[2821] = v[1900] + v[1901];
	v[2814] = v[1912] + v[2821];
	v[293] = v[281] * v[723];
	v[1916] = v[293] * v[4740];
	v[295] = v[281] * v[722];
	v[5249] = v[295] / v[442];
	v[1896] = v[295] * v[4742];
	v[296] = 1e0 + v[4741] * v[719];
	v[1897] = v[296] * v[4743];
	v[2819] = v[1896] + v[1897] + v[1916];
	v[2816] = -v[1916] + v[2819];
	v[2810] = -v[1896] + v[2819];
	v[297] = -(v[114] * v[4741]);
	v[4991] = 2e0 * v[297];
	v[448] = -(v[281] * v[297]);
	v[298] = v[113] * v[4741];
	v[4992] = 2e0 * v[298];
	v[460] = -(v[297] * v[298]);
	v[457] = -(v[281] * v[298]);
	v[299] = -(v[112] * v[4741]);
	v[4993] = 2e0 * v[299];
	v[461] = -(v[281] * v[299]);
	v[4768] = v[460] - v[461];
	v[458] = v[297] * v[299];
	v[4766] = -v[457] + v[458];
	v[449] = -(v[298] * v[299]);
	v[4764] = -v[448] + v[449];
	v[300] = v[148] * v[265] + v[151] * v[266] + v[154] * v[267];
	v[301] = v[149] * v[265] + v[152] * v[266] + v[155] * v[267];
	v[303] = v[148] * v[269] + v[151] * v[271] + v[154] * v[272];
	v[304] = v[149] * v[269] + v[152] * v[271] + v[155] * v[272];
	v[306] = v[148] * v[274] + v[151] * v[276] + v[154] * v[277];
	v[307] = v[149] * v[274] + v[152] * v[276] + v[155] * v[277];
	v[309] = v[157] * v[284] + v[160] * v[285] + v[163] * v[286];
	v[310] = v[158] * v[284] + v[161] * v[285] + v[164] * v[286];
	v[312] = v[157] * v[288] + v[160] * v[290] + v[163] * v[291];
	v[313] = v[158] * v[288] + v[161] * v[290] + v[164] * v[291];
	v[315] = v[157] * v[293] + v[160] * v[295] + v[163] * v[296];
	v[316] = v[158] * v[293] + v[161] * v[295] + v[164] * v[296];
	v[318] = v[103] + v[142];
	v[319] = v[104] + v[143];
	v[320] = v[105] + v[144];
	v[321] = v[109] + v[145];
	v[322] = v[110] + v[146];
	v[323] = v[111] + v[147];
	v[324] = (1e0 - v[7]) / 2e0;
	v[325] = (1e0 + v[7]) / 2e0;
	v[326] = (1e0 - v[3]) / 2e0;
	v[327] = (1e0 + v[3]) / 2e0;
	v[9988] = 0e0;
	v[9989] = 0e0;
	v[9990] = v[239];
	v[9991] = 0e0;
	v[9992] = 0e0;
	v[9993] = 0e0;
	v[9994] = 0e0;
	v[9995] = 0e0;
	v[9996] = v[240];
	v[9997] = 0e0;
	v[9998] = 0e0;
	v[9999] = 0e0;
	v[10000] = 0e0;
	v[10001] = 0e0;
	v[10002] = -v[326];
	v[10003] = 0e0;
	v[10004] = 0e0;
	v[10005] = 0e0;
	v[10006] = 0e0;
	v[10007] = 0e0;
	v[10008] = -v[327];
	v[10009] = 0e0;
	v[10010] = 0e0;
	v[10011] = 0e0;
	v[10012] = 0e0;
	v[10013] = v[239];
	v[10014] = 0e0;
	v[10015] = 0e0;
	v[10016] = 0e0;
	v[10017] = 0e0;
	v[10018] = 0e0;
	v[10019] = v[240];
	v[10020] = 0e0;
	v[10021] = 0e0;
	v[10022] = 0e0;
	v[10023] = 0e0;
	v[10024] = 0e0;
	v[10025] = -v[326];
	v[10026] = 0e0;
	v[10027] = 0e0;
	v[10028] = 0e0;
	v[10029] = 0e0;
	v[10030] = 0e0;
	v[10031] = -v[327];
	v[10032] = 0e0;
	v[10033] = 0e0;
	v[10034] = 0e0;
	v[10035] = 0e0;
	v[10036] = v[239];
	v[10037] = 0e0;
	v[10038] = 0e0;
	v[10039] = 0e0;
	v[10040] = 0e0;
	v[10041] = 0e0;
	v[10042] = v[240];
	v[10043] = 0e0;
	v[10044] = 0e0;
	v[10045] = 0e0;
	v[10046] = 0e0;
	v[10047] = 0e0;
	v[10048] = -v[326];
	v[10049] = 0e0;
	v[10050] = 0e0;
	v[10051] = 0e0;
	v[10052] = 0e0;
	v[10053] = 0e0;
	v[10054] = -v[327];
	v[10055] = 0e0;
	v[10056] = 0e0;
	v[10057] = 0e0;
	v[10058] = 0e0;
	v[10059] = 0e0;
	v[4034] = v[306] * v[326] + v[315] * v[327];
	v[4033] = v[303] * v[326] + v[312] * v[327];
	v[4032] = v[300] * v[326] + v[309] * v[327];
	v[4030] = v[307] * v[326] + v[316] * v[327];
	v[4029] = v[304] * v[326] + v[313] * v[327];
	v[4028] = v[301] * v[326] + v[310] * v[327];
	v[328] = 2e0 / (*eB);
	v[4744] = -1e0 + v[328];
	v[3968] = Power(v[532], v[4744]);
	v[5223] = -(v[1465] * v[3968]);
	v[3967] = -2e0 + v[328];
	v[3966] = Power(v[534], v[4744]);
	v[5226] = v[1463] * v[3966];
	v[1468] = v[328] * (v[1464] * v[3966] + v[1466] * v[3968]);
	v[4699] = v[1466] * v[4744] * Power(v[532], v[3967]);
	v[4697] = v[1464] * v[4744] * Power(v[534], v[3967]);
	v[531] = Power(v[532], v[328]) + Power(v[534], v[328]);
	v[3624] = -1e0 - 1e0 / v[328];
	v[5265] = -((v[1468] * v[3624] * Power(v[531], -1e0 + v[3624])) / v[328]);
	v[3623] = Power(v[531], v[3624]);
	v[5088] = -(v[3623] / v[328]);
	v[535] = v[1468] * v[5088];
	v[329] = 1e0 / Power(Power(fabs(v[330]), v[328]) + Power(fabs(v[332]), v[328]), 1 / v[328]);
	v[331] = v[139] * v[329] * v[330];
	v[333] = v[140] * v[329] * v[332];
	v[335] = 1e0 / Power(v[531], 1 / v[328]);
	v[5224] = v[140] * v[335];
	v[4747] = v[139] * v[335];
	v[4262] = v[1565] * v[5224];
	v[4260] = -(v[1462] * v[4747]);
	v[537] = v[4262] + v[4748] * v[535];
	v[536] = v[4260] + v[5087] * v[535];
	v[543] = v[300] * v[536] + v[301] * v[537];
	v[542] = v[309] * v[536] + v[310] * v[537];
	v[544] = v[327] * v[542] + v[326] * v[543];
	v[541] = v[303] * v[536] + v[304] * v[537];
	v[540] = v[312] * v[536] + v[313] * v[537];
	v[545] = v[327] * v[540] + v[326] * v[541];
	v[539] = v[306] * v[536] + v[307] * v[537];
	v[538] = v[315] * v[536] + v[316] * v[537];
	v[546] = v[327] * v[538] + v[326] * v[539];
	v[337] = v[1565] * v[4747];
	v[339] = v[335] * v[4748];
	v[3631] = v[163] * v[337] + v[164] * v[339];
	v[3630] = v[160] * v[337] + v[161] * v[339];
	v[3629] = v[157] * v[337] + v[158] * v[339];
	v[3628] = v[154] * v[337] + v[155] * v[339];
	v[3627] = v[151] * v[337] + v[152] * v[339];
	v[3626] = v[148] * v[337] + v[149] * v[339];
	v[527] = v[323] + v[315] * v[337] + v[316] * v[339];
	v[526] = v[320] + v[306] * v[337] + v[307] * v[339];
	v[528] = (-v[526] + v[527]) / 2e0;
	v[524] = v[322] + v[312] * v[337] + v[313] * v[339];
	v[523] = v[319] + v[303] * v[337] + v[304] * v[339];
	v[525] = (-v[523] + v[524]) / 2e0;
	v[521] = v[321] + v[309] * v[337] + v[310] * v[339];
	v[520] = v[318] + v[300] * v[337] + v[301] * v[339];
	v[522] = (-v[520] + v[521]) / 2e0;
	v[349] = duiA[0] * v[20] + dduiA[0] * v[21] + v[19] * v[40];
	v[350] = duiA[1] * v[20] + dduiA[1] * v[21] + v[19] * v[41];
	v[351] = duiA[2] * v[20] + dduiA[2] * v[21] + v[19] * v[42];
	v[4862] = v[239] * v[351];
	v[352] = duiA[6] * v[20] + dduiA[6] * v[21] + v[19] * v[46];
	v[353] = duiA[7] * v[20] + dduiA[7] * v[21] + v[19] * v[47];
	v[354] = duiA[8] * v[20] + dduiA[8] * v[21] + v[19] * v[48];
	v[4863] = v[240] * v[354];
	v[355] = v[103] * v[19] + duiB[0] * v[20] + dduiB[0] * v[21];
	v[356] = v[104] * v[19] + duiB[1] * v[20] + dduiB[1] * v[21];
	v[357] = v[105] * v[19] + duiB[2] * v[20] + dduiB[2] * v[21];
	v[4864] = -(v[326] * v[357]);
	v[358] = v[109] * v[19] + duiB[6] * v[20] + dduiB[6] * v[21];
	v[359] = v[110] * v[19] + duiB[7] * v[20] + dduiB[7] * v[21];
	v[360] = v[111] * v[19] + duiB[8] * v[20] + dduiB[8] * v[21];
	v[4865] = -(v[327] * v[360]);
	v[362] = v[379] + v[380];
	v[4752] = v[362] / v[364];
	v[363] = v[370] + v[371];
	v[4750] = v[363] / v[364];
	v[365] = (v[193] * v[193]) + v[364];
	v[5158] = v[190] * v[362] + v[180] * v[365];
	v[5157] = v[184] * v[363] + v[179] * v[365];
	v[3976] = v[365] / v[364];
	v[2212] = v[4751] / v[364];
	v[2210] = v[4749] / v[364];
	v[8006] = 0e0;
	v[8007] = 0e0;
	v[8008] = 0e0;
	v[8009] = 0e0;
	v[8010] = v[19] * v[2210];
	v[8011] = v[19] * v[2212];
	v[8012] = 0e0;
	v[8013] = 0e0;
	v[8014] = 0e0;
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
	v[8027] = 0e0;
	v[8028] = 0e0;
	v[8029] = 0e0;
	v[366] = v[2210] * v[373] + v[2212] * v[375] + v[367] * v[3976];
	v[368] = (v[192] * v[192]) + v[364];
	v[5155] = v[182] * v[368] + v[178] * v[4749];
	v[374] = v[382] + v[383];
	v[5156] = v[185] * v[368] + v[190] * v[374];
	v[3978] = v[368] / v[364];
	v[2211] = v[4753] / v[364];
	v[8054] = 0e0;
	v[8055] = 0e0;
	v[8056] = 0e0;
	v[8057] = v[19] * v[4750];
	v[8058] = 0e0;
	v[8059] = v[19] * v[2211];
	v[8060] = 0e0;
	v[8061] = 0e0;
	v[8062] = 0e0;
	v[8063] = 0e0;
	v[8064] = 0e0;
	v[8065] = 0e0;
	v[8066] = 0e0;
	v[8067] = 0e0;
	v[8068] = 0e0;
	v[8069] = 0e0;
	v[8070] = 0e0;
	v[8071] = 0e0;
	v[8072] = 0e0;
	v[8073] = 0e0;
	v[8074] = 0e0;
	v[8075] = 0e0;
	v[8076] = 0e0;
	v[8077] = 0e0;
	v[376] = v[2211] * v[375] + v[373] * v[3978] + v[367] * v[4750];
	v[377] = (v[191] * v[191]) + v[364];
	v[5153] = v[187] * v[377] + v[178] * v[4751];
	v[5154] = v[189] * v[377] + v[184] * v[4753];
	v[3980] = v[377] / v[364];
	v[2209] = v[374] / v[364];
	v[8102] = 0e0;
	v[8103] = 0e0;
	v[8104] = 0e0;
	v[8105] = v[19] * v[4752];
	v[8106] = v[19] * v[2209];
	v[8107] = 0e0;
	v[8108] = 0e0;
	v[8109] = 0e0;
	v[8110] = 0e0;
	v[8111] = 0e0;
	v[8112] = 0e0;
	v[8113] = 0e0;
	v[8114] = 0e0;
	v[8115] = 0e0;
	v[8116] = 0e0;
	v[8117] = 0e0;
	v[8118] = 0e0;
	v[8119] = 0e0;
	v[8120] = 0e0;
	v[8121] = 0e0;
	v[8122] = 0e0;
	v[8123] = 0e0;
	v[8124] = 0e0;
	v[8125] = 0e0;
	v[386] = v[2209] * v[373] + v[375] * v[3980] + v[367] * v[4752];
	v[388] = v[405] + v[406];
	v[4757] = v[388] / v[390];
	v[389] = v[396] + v[397];
	v[4755] = v[389] / v[390];
	v[391] = (v[212] * v[212]) + v[390];
	v[5145] = v[209] * v[388] + v[199] * v[391];
	v[5144] = v[203] * v[389] + v[198] * v[391];
	v[3982] = v[391] / v[390];
	v[2218] = v[4756] / v[390];
	v[2216] = v[4754] / v[390];
	v[8150] = 0e0;
	v[8151] = 0e0;
	v[8152] = 0e0;
	v[8153] = 0e0;
	v[8154] = 0e0;
	v[8155] = 0e0;
	v[8156] = 0e0;
	v[8157] = 0e0;
	v[8158] = 0e0;
	v[8159] = 0e0;
	v[8160] = v[19] * v[2216];
	v[8161] = v[19] * v[2218];
	v[8162] = 0e0;
	v[8163] = 0e0;
	v[8164] = 0e0;
	v[8165] = 0e0;
	v[8166] = 0e0;
	v[8167] = 0e0;
	v[8168] = 0e0;
	v[8169] = 0e0;
	v[8170] = 0e0;
	v[8171] = 0e0;
	v[8172] = 0e0;
	v[8173] = 0e0;
	v[392] = v[393] * v[3982] + v[2216] * v[399] + v[2218] * v[401];
	v[394] = (v[211] * v[211]) + v[390];
	v[5142] = v[201] * v[394] + v[197] * v[4754];
	v[400] = v[408] + v[409];
	v[5143] = v[204] * v[394] + v[209] * v[400];
	v[3984] = v[394] / v[390];
	v[2217] = v[4758] / v[390];
	v[8198] = 0e0;
	v[8199] = 0e0;
	v[8200] = 0e0;
	v[8201] = 0e0;
	v[8202] = 0e0;
	v[8203] = 0e0;
	v[8204] = 0e0;
	v[8205] = 0e0;
	v[8206] = 0e0;
	v[8207] = v[19] * v[4755];
	v[8208] = 0e0;
	v[8209] = v[19] * v[2217];
	v[8210] = 0e0;
	v[8211] = 0e0;
	v[8212] = 0e0;
	v[8213] = 0e0;
	v[8214] = 0e0;
	v[8215] = 0e0;
	v[8216] = 0e0;
	v[8217] = 0e0;
	v[8218] = 0e0;
	v[8219] = 0e0;
	v[8220] = 0e0;
	v[8221] = 0e0;
	v[402] = v[3984] * v[399] + v[2217] * v[401] + v[393] * v[4755];
	v[403] = (v[210] * v[210]) + v[390];
	v[5140] = v[206] * v[403] + v[197] * v[4756];
	v[5141] = v[208] * v[403] + v[203] * v[4758];
	v[3986] = v[403] / v[390];
	v[2215] = v[400] / v[390];
	v[8246] = 0e0;
	v[8247] = 0e0;
	v[8248] = 0e0;
	v[8249] = 0e0;
	v[8250] = 0e0;
	v[8251] = 0e0;
	v[8252] = 0e0;
	v[8253] = 0e0;
	v[8254] = 0e0;
	v[8255] = v[19] * v[4757];
	v[8256] = v[19] * v[2215];
	v[8257] = 0e0;
	v[8258] = 0e0;
	v[8259] = 0e0;
	v[8260] = 0e0;
	v[8261] = 0e0;
	v[8262] = 0e0;
	v[8263] = 0e0;
	v[8264] = 0e0;
	v[8265] = 0e0;
	v[8266] = 0e0;
	v[8267] = 0e0;
	v[8268] = 0e0;
	v[8269] = 0e0;
	v[412] = v[2215] * v[399] + v[3986] * v[401] + v[393] * v[4757];
	v[414] = v[431] + v[432];
	v[4762] = v[414] / v[416];
	v[415] = v[422] + v[423];
	v[4760] = v[415] / v[416];
	v[417] = (v[280] * v[280]) + v[416];
	v[5132] = v[277] * v[414] + v[267] * v[417];
	v[5131] = v[271] * v[415] + v[266] * v[417];
	v[3988] = v[417] / v[416];
	v[2224] = v[4761] / v[416];
	v[2222] = v[4759] / v[416];
	v[8294] = 0e0;
	v[8295] = 0e0;
	v[8296] = 0e0;
	v[8297] = 0e0;
	v[8298] = 0e0;
	v[8299] = 0e0;
	v[8300] = 0e0;
	v[8301] = 0e0;
	v[8302] = 0e0;
	v[8303] = 0e0;
	v[8304] = 0e0;
	v[8305] = 0e0;
	v[8306] = 0e0;
	v[8307] = 0e0;
	v[8308] = 0e0;
	v[8309] = 0e0;
	v[8310] = v[19] * v[2222];
	v[8311] = v[19] * v[2224];
	v[8312] = 0e0;
	v[8313] = 0e0;
	v[8314] = 0e0;
	v[8315] = 0e0;
	v[8316] = 0e0;
	v[8317] = 0e0;
	v[418] = v[3988] * v[419] + v[2222] * v[425] + v[2224] * v[427];
	v[420] = (v[279] * v[279]) + v[416];
	v[5129] = v[269] * v[420] + v[265] * v[4759];
	v[426] = v[434] + v[435];
	v[5130] = v[272] * v[420] + v[277] * v[426];
	v[3990] = v[420] / v[416];
	v[2223] = v[4763] / v[416];
	v[8342] = 0e0;
	v[8343] = 0e0;
	v[8344] = 0e0;
	v[8345] = 0e0;
	v[8346] = 0e0;
	v[8347] = 0e0;
	v[8348] = 0e0;
	v[8349] = 0e0;
	v[8350] = 0e0;
	v[8351] = 0e0;
	v[8352] = 0e0;
	v[8353] = 0e0;
	v[8354] = 0e0;
	v[8355] = 0e0;
	v[8356] = 0e0;
	v[8357] = v[19] * v[4760];
	v[8358] = 0e0;
	v[8359] = v[19] * v[2223];
	v[8360] = 0e0;
	v[8361] = 0e0;
	v[8362] = 0e0;
	v[8363] = 0e0;
	v[8364] = 0e0;
	v[8365] = 0e0;
	v[428] = v[3990] * v[425] + v[2223] * v[427] + v[419] * v[4760];
	v[429] = (v[278] * v[278]) + v[416];
	v[5127] = v[274] * v[429] + v[265] * v[4761];
	v[5128] = v[276] * v[429] + v[271] * v[4763];
	v[3992] = v[429] / v[416];
	v[2221] = v[426] / v[416];
	v[8390] = 0e0;
	v[8391] = 0e0;
	v[8392] = 0e0;
	v[8393] = 0e0;
	v[8394] = 0e0;
	v[8395] = 0e0;
	v[8396] = 0e0;
	v[8397] = 0e0;
	v[8398] = 0e0;
	v[8399] = 0e0;
	v[8400] = 0e0;
	v[8401] = 0e0;
	v[8402] = 0e0;
	v[8403] = 0e0;
	v[8404] = 0e0;
	v[8405] = v[19] * v[4762];
	v[8406] = v[19] * v[2221];
	v[8407] = 0e0;
	v[8408] = 0e0;
	v[8409] = 0e0;
	v[8410] = 0e0;
	v[8411] = 0e0;
	v[8412] = 0e0;
	v[8413] = 0e0;
	v[438] = v[2221] * v[425] + v[3992] * v[427] + v[419] * v[4762];
	v[440] = v[457] + v[458];
	v[4767] = v[440] / v[442];
	v[441] = v[448] + v[449];
	v[4765] = v[441] / v[442];
	v[443] = (v[299] * v[299]) + v[442];
	v[5119] = v[296] * v[440] + v[286] * v[443];
	v[5118] = v[290] * v[441] + v[285] * v[443];
	v[3994] = v[443] / v[442];
	v[2230] = v[4766] / v[442];
	v[2228] = v[4764] / v[442];
	v[8438] = 0e0;
	v[8439] = 0e0;
	v[8440] = 0e0;
	v[8441] = 0e0;
	v[8442] = 0e0;
	v[8443] = 0e0;
	v[8444] = 0e0;
	v[8445] = 0e0;
	v[8446] = 0e0;
	v[8447] = 0e0;
	v[8448] = 0e0;
	v[8449] = 0e0;
	v[8450] = 0e0;
	v[8451] = 0e0;
	v[8452] = 0e0;
	v[8453] = 0e0;
	v[8454] = 0e0;
	v[8455] = 0e0;
	v[8456] = 0e0;
	v[8457] = 0e0;
	v[8458] = 0e0;
	v[8459] = 0e0;
	v[8460] = v[19] * v[2228];
	v[8461] = v[19] * v[2230];
	v[444] = v[3994] * v[445] + v[2228] * v[451] + v[2230] * v[453];
	v[446] = (v[298] * v[298]) + v[442];
	v[5116] = v[288] * v[446] + v[284] * v[4764];
	v[452] = v[460] + v[461];
	v[5117] = v[291] * v[446] + v[296] * v[452];
	v[3996] = v[446] / v[442];
	v[2229] = v[4768] / v[442];
	v[8486] = 0e0;
	v[8487] = 0e0;
	v[8488] = 0e0;
	v[8489] = 0e0;
	v[8490] = 0e0;
	v[8491] = 0e0;
	v[8492] = 0e0;
	v[8493] = 0e0;
	v[8494] = 0e0;
	v[8495] = 0e0;
	v[8496] = 0e0;
	v[8497] = 0e0;
	v[8498] = 0e0;
	v[8499] = 0e0;
	v[8500] = 0e0;
	v[8501] = 0e0;
	v[8502] = 0e0;
	v[8503] = 0e0;
	v[8504] = 0e0;
	v[8505] = 0e0;
	v[8506] = 0e0;
	v[8507] = v[19] * v[4765];
	v[8508] = 0e0;
	v[8509] = v[19] * v[2229];
	v[454] = v[3996] * v[451] + v[2229] * v[453] + v[445] * v[4765];
	v[455] = (v[297] * v[297]) + v[442];
	v[5114] = v[293] * v[455] + v[284] * v[4766];
	v[5115] = v[295] * v[455] + v[290] * v[4768];
	v[3998] = v[455] / v[442];
	v[2227] = v[452] / v[442];
	v[8534] = 0e0;
	v[8535] = 0e0;
	v[8536] = 0e0;
	v[8537] = 0e0;
	v[8538] = 0e0;
	v[8539] = 0e0;
	v[8540] = 0e0;
	v[8541] = 0e0;
	v[8542] = 0e0;
	v[8543] = 0e0;
	v[8544] = 0e0;
	v[8545] = 0e0;
	v[8546] = 0e0;
	v[8547] = 0e0;
	v[8548] = 0e0;
	v[8549] = 0e0;
	v[8550] = 0e0;
	v[8551] = 0e0;
	v[8552] = 0e0;
	v[8553] = 0e0;
	v[8554] = 0e0;
	v[8555] = v[19] * v[4767];
	v[8556] = v[19] * v[2227];
	v[8557] = 0e0;
	v[464] = v[2227] * v[451] + v[3998] * v[453] + v[445] * v[4767];
	v[465] = v[239] * v[492] + v[240] * v[493] - v[326] * v[520] - v[327] * v[521];
	v[4791] = v[240] * v[465];
	v[4790] = v[239] * v[465];
	v[4779] = -(v[327] * v[465]);
	v[4778] = -(v[326] * v[465]);
	v[588] = v[465] / 2e0;
	v[466] = v[239] * v[495] + v[240] * v[496] - v[326] * v[523] - v[327] * v[524];
	v[4793] = v[240] * v[466];
	v[4792] = v[239] * v[466];
	v[4781] = -(v[327] * v[466]);
	v[4780] = -(v[326] * v[466]);
	v[571] = v[466] / 2e0;
	v[467] = v[239] * v[498] + v[240] * v[499] - v[326] * v[526] - v[327] * v[527];
	v[4795] = v[240] * v[467];
	v[4794] = v[239] * v[467];
	v[4783] = -(v[327] * v[467]);
	v[4782] = -(v[326] * v[467]);
	v[1475] = sqrt((v[465] * v[465]) + (v[466] * v[466]) + (v[467] * v[467]));
	v[5245] = v[12] * v[13] * Power(v[1475], -1e0 + v[13]);
	v[3469] = 1e0 / (v[1475] * v[1475]);
	v[3165] = v[12] * Power(v[1475], v[13]);
	v[1347] = -(v[327] * v[3631]);
	v[5274] = v[1347] / 2e0;
	v[1348] = -(v[327] * v[3630]);
	v[5275] = v[1348] / 2e0;
	v[1349] = -(v[327] * v[3629]);
	v[5276] = v[1349] / 2e0;
	v[1350] = -(v[326] * v[3628]);
	v[5277] = v[1350] / 2e0;
	v[1351] = -(v[326] * v[3627]);
	v[5278] = v[1351] / 2e0;
	v[1352] = -(v[326] * v[3626]);
	v[5279] = v[1352] / 2e0;
	v[1353] = v[240] * v[3765];
	v[5293] = v[1353] / 2e0;
	v[1354] = v[240] * v[3764];
	v[5294] = v[1354] / 2e0;
	v[1355] = v[240] * v[3763];
	v[5295] = v[1355] / 2e0;
	v[1356] = v[239] * v[3762];
	v[5296] = v[1356] / 2e0;
	v[1357] = v[239] * v[3761];
	v[5297] = v[1357] / 2e0;
	v[1358] = v[239] * v[3760];
	v[5298] = v[1358] / 2e0;
	v[1359] = v[1347] * v[281];
	v[1360] = v[1348] * v[281];
	v[1361] = v[1349] * v[281];
	v[1362] = v[1349] * v[4786] + v[1347] * v[732] + v[1348] * v[737];
	v[4825] = v[1362] * v[749];
	v[1363] = v[1348] * v[4785] + v[1347] * v[724] + v[1349] * v[731];
	v[4826] = v[1363] * v[749];
	v[1364] = v[1347] * v[4784] + v[1348] * v[722] + v[1349] * v[723];
	v[4827] = v[1364] * v[749];
	v[1365] = v[1350] * v[262];
	v[1366] = v[1351] * v[262];
	v[1367] = v[1352] * v[262];
	v[1368] = v[1352] * v[4789] + v[1350] * v[763] + v[1351] * v[768];
	v[4831] = v[1368] * v[780];
	v[1369] = v[1351] * v[4788] + v[1350] * v[755] + v[1352] * v[762];
	v[4832] = v[1369] * v[780];
	v[1370] = v[1350] * v[4787] + v[1351] * v[753] + v[1352] * v[754];
	v[4833] = v[1370] * v[780];
	v[1371] = v[1353] * v[194];
	v[1372] = v[1354] * v[194];
	v[1373] = v[1355] * v[194];
	v[1374] = v[1355] * v[4798] + v[1353] * v[914] + v[1354] * v[919];
	v[4837] = v[1374] * v[931];
	v[1375] = v[1354] * v[4797] + v[1353] * v[906] + v[1355] * v[913];
	v[4838] = v[1375] * v[931];
	v[1376] = v[1353] * v[4796] + v[1354] * v[904] + v[1355] * v[905];
	v[4839] = v[1376] * v[931];
	v[1377] = v[1356] * v[175];
	v[1378] = v[1357] * v[175];
	v[1379] = v[1358] * v[175];
	v[1380] = v[1358] * v[4801] + v[1356] * v[945] + v[1357] * v[950];
	v[4843] = v[1380] * v[962];
	v[1381] = v[1357] * v[4800] + v[1356] * v[937] + v[1358] * v[944];
	v[4846] = v[1381] * v[962];
	v[1382] = v[1356] * v[4799] + v[1357] * v[935] + v[1358] * v[936];
	v[4859] = v[1382] * v[962];
	v[554] = v[467] / 2e0;
	v[468] = -(v[324] * (v[142] + v[148] * v[331] + v[149] * v[333])) - v[325] * (v[145] + v[157] * v[331] + v[158] * v[333])
		+ v[237] * (v[79] + v[244] * v[85] + v[246] * v[86]) + v[238] * (v[82] + v[244] * v[94] + v[246] * v[95]);
	v[469] = -(v[324] * (v[143] + v[151] * v[331] + v[152] * v[333])) - v[325] * (v[146] + v[160] * v[331] + v[161] * v[333])
		+ v[237] * (v[80] + v[244] * v[88] + v[246] * v[89]) + v[238] * (v[83] + v[244] * v[97] + v[246] * v[98]);
	v[470] = -(v[324] * (v[144] + v[154] * v[331] + v[155] * v[333])) - v[325] * (v[147] + v[163] * v[331] + v[164] * v[333])
		+ v[238] * (v[100] * v[244] + v[101] * v[246] + v[84]) + v[237] * (v[81] + v[244] * v[91] + v[246] * v[92]);
	if (v[1475] > 0.1e-7) { v01 = 1e0 / v[1475]; v02 = (-(v01 / v[1475])); v03 = (2e0 * v01) / (v[1475] * v[1475]); }
	else {
		v01 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[1475]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[1475]) *
			(0.2399999997e10 - 0.1199999994e18 * v[1475] - 0.3e17 * (v[1475] * v[1475]))));
		v02 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[1475] + 0.6e25 * Power(v[1475], 3)
			+ 0.1799999982e26 * (v[1475] * v[1475]));
		v03 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[1475] - 0.3e17 * (v[1475] * v[1475]));
	};
	v[477] = v03;
	v[478] = v02;
	v[479] = v01;
	v[480] = v[465] * v[479];
	v[4777] = -(v[355] * v[480]);
	v[4776] = -(v[358] * v[480]);
	v[4775] = v[349] * v[480];
	v[4774] = v[352] * v[480];
	v[1383] = (v[480] * v[480]);
	v[481] = v[466] * v[479];
	v[5090] = v[480] * v[481];
	v[4772] = v[239] * v[481];
	v[4771] = v[240] * v[481];
	v[4770] = -(v[326] * v[481]);
	v[4769] = -(v[327] * v[481]);
	v[4275] = v[17] * v[5090];
	v[4893] = v[359] * v[4769] + v[356] * v[4770] + v[353] * v[4771] + v[350] * v[4772] + v[240] * v[4774] + v[239] * v[4775]
		+ v[327] * v[4776] + v[326] * v[4777];
	v[1413] = (v[481] * v[481]);
	v[1412] = v[4769] * v[480];
	v[1411] = v[4770] * v[480];
	v[1410] = v[4771] * v[480];
	v[1409] = v[4772] * v[480];
	v[482] = v[467] * v[479];
	v[4773] = v[17] * v[482];
	v[4315] = v[464] * v[4773];
	v[4311] = v[444] * v[4773];
	v[4309] = v[454] * v[4773];
	v[4303] = v[438] * v[4773];
	v[4299] = v[418] * v[4773];
	v[4297] = v[428] * v[4773];
	v[4291] = v[412] * v[4773];
	v[4287] = v[392] * v[4773];
	v[4285] = v[402] * v[4773];
	v[4277] = v[386] * v[4773];
	v[4271] = v[366] * v[4773];
	v[4269] = v[376] * v[4773];
	v[1455] = (v[482] * v[482]);
	v[4608] = v[1455] * v[354] + (v[4774] + v[353] * v[481]) * v[482];
	v[4607] = v[1455] * v[351] + (v[4775] + v[350] * v[481]) * v[482];
	v[4527] = -(v[1455] * v[360]) + (v[4776] - v[359] * v[481]) * v[482];
	v[4526] = -(v[1455] * v[357]) + (v[4777] - v[356] * v[481]) * v[482];
	v[3251] = v[482] * (v[4862] + v[4863] + v[4864] + v[4865]);
	v[483] = sqrt((v[468] * v[468]) + (v[469] * v[469]) + (v[470] * v[470]));
	if (v[483] > 0.1e-7) { v04 = 1e0 / v[483]; v05 = (-(v04 / v[483])); v06 = (2e0 * v04) / (v[483] * v[483]); }
	else {
		v04 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[483]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[483]) * (0.2399999997e10
			- 0.1199999994e18 * v[483] - 0.3e17 * (v[483] * v[483]))));
		v05 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[483] + 0.6e25 * Power(v[483], 3)
			+ 0.1799999982e26 * (v[483] * v[483]));
		v06 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[483] - 0.3e17 * (v[483] * v[483]));
	};
	v[488] = v04;
	v[489] = v[468] * v[488];
	v[490] = v[469] * v[488];
	v[491] = v[470] * v[488];
	v[548] = -(v[326] * v[500]);
	v[549] = -(v[326] * v[518]);
	v[550] = v[326] * v[528] + v[554];
	v[551] = v[326] * v[546];
	v[552] = -(v[327] * v[500]);
	v[553] = -(v[327] * v[518]);
	v[555] = v[327] * v[528] - v[554];
	v[556] = v[327] * v[546];
	v[557] = v[239] * v[500] - v[554];
	v[558] = v[239] * v[518];
	v[559] = -(v[239] * v[528]);
	v[560] = -(v[239] * v[546]);
	v[561] = v[240] * v[500] + v[554];
	v[562] = v[240] * v[518];
	v[563] = -(v[240] * v[528]);
	v[564] = -(v[240] * v[546]);
	v[565] = -(v[326] * v[497]);
	v[566] = -(v[326] * v[517]);
	v[567] = v[326] * v[525] + v[571];
	v[568] = v[326] * v[545];
	v[569] = -(v[327] * v[497]);
	v[570] = -(v[327] * v[517]);
	v[572] = v[327] * v[525] - v[571];
	v[573] = v[327] * v[545];
	v[574] = v[239] * v[497] - v[571];
	v[575] = v[239] * v[517];
	v[576] = -(v[239] * v[525]);
	v[577] = -(v[239] * v[545]);
	v[578] = v[240] * v[497] + v[571];
	v[579] = v[240] * v[517];
	v[580] = -(v[240] * v[525]);
	v[581] = -(v[240] * v[545]);
	v[582] = -(v[326] * v[494]);
	v[583] = -(v[326] * v[516]);
	v[584] = v[326] * v[522] + v[588];
	v[585] = v[326] * v[544];
	v[586] = -(v[327] * v[494]);
	v[587] = -(v[327] * v[516]);
	v[589] = v[327] * v[522] - v[588];
	v[590] = v[327] * v[544];
	v[591] = v[239] * v[494] - v[588];
	v[592] = v[239] * v[516];
	v[593] = -(v[239] * v[522]);
	v[594] = -(v[239] * v[544]);
	v[595] = v[240] * v[494] + v[588];
	v[596] = v[240] * v[516];
	v[597] = -(v[240] * v[522]);
	v[598] = -(v[240] * v[544]);
	v[599] = v[339] * v[582];
	v[600] = v[339] * v[583];
	v[601] = v[339] * v[584];
	v[602] = v[4778] * v[537] + v[339] * v[585];
	v[603] = v[337] * v[582];
	v[604] = v[337] * v[583];
	v[605] = v[337] * v[584];
	v[606] = v[4778] * v[536] + v[337] * v[585];
	v[607] = v[339] * v[586];
	v[608] = v[339] * v[587];
	v[609] = v[339] * v[589];
	v[610] = v[4779] * v[537] + v[339] * v[590];
	v[611] = v[337] * v[586];
	v[612] = v[337] * v[587];
	v[613] = v[337] * v[589];
	v[614] = v[4779] * v[536] + v[337] * v[590];
	v[615] = v[339] * v[565];
	v[616] = v[339] * v[566];
	v[617] = v[339] * v[567];
	v[618] = v[4780] * v[537] + v[339] * v[568];
	v[619] = v[337] * v[565];
	v[620] = v[337] * v[566];
	v[621] = v[337] * v[567];
	v[622] = v[4780] * v[536] + v[337] * v[568];
	v[623] = v[339] * v[569];
	v[624] = v[339] * v[570];
	v[625] = v[339] * v[572];
	v[626] = v[4781] * v[537] + v[339] * v[573];
	v[627] = v[337] * v[569];
	v[628] = v[337] * v[570];
	v[629] = v[337] * v[572];
	v[630] = v[4781] * v[536] + v[337] * v[573];
	v[631] = v[339] * v[548];
	v[632] = v[339] * v[549];
	v[633] = v[339] * v[550];
	v[634] = v[4782] * v[537] + v[339] * v[551];
	v[635] = v[337] * v[548];
	v[636] = v[337] * v[549];
	v[637] = v[337] * v[550];
	v[638] = v[4782] * v[536] + v[337] * v[551];
	v[639] = v[339] * v[552];
	v[640] = v[339] * v[553];
	v[641] = v[339] * v[555];
	v[642] = v[4783] * v[537] + v[339] * v[556];
	v[643] = v[337] * v[552];
	v[644] = v[337] * v[553];
	v[645] = v[337] * v[555];
	v[646] = v[4783] * v[536] + v[337] * v[556];
	v[647] = v[164] * v[639] + v[163] * v[643];
	v[648] = v[164] * v[640] + v[163] * v[644];
	v[649] = v[164] * v[641] + v[163] * v[645];
	v[650] = v[164] * v[642] + v[163] * v[646];
	v[651] = v[161] * v[639] + v[160] * v[643];
	v[1026] = v[281] * v[651];
	v[652] = v[161] * v[640] + v[160] * v[644];
	v[1078] = v[281] * v[652];
	v[653] = v[161] * v[641] + v[160] * v[645];
	v[1130] = v[281] * v[653];
	v[654] = v[161] * v[642] + v[160] * v[646];
	v[1182] = v[281] * v[654];
	v[655] = v[158] * v[639] + v[157] * v[643];
	v[1029] = v[281] * v[655];
	v[656] = v[158] * v[640] + v[157] * v[644];
	v[1081] = v[281] * v[656];
	v[657] = v[158] * v[641] + v[157] * v[645];
	v[1133] = v[281] * v[657];
	v[658] = v[158] * v[642] + v[157] * v[646];
	v[1185] = v[281] * v[658];
	v[659] = v[164] * v[623] + v[163] * v[627];
	v[1027] = v[281] * v[659];
	v[660] = v[164] * v[624] + v[163] * v[628];
	v[1079] = v[281] * v[660];
	v[661] = v[164] * v[625] + v[163] * v[629];
	v[1131] = v[281] * v[661];
	v[662] = v[164] * v[626] + v[163] * v[630];
	v[1183] = v[281] * v[662];
	v[663] = v[161] * v[623] + v[160] * v[627];
	v[664] = v[161] * v[624] + v[160] * v[628];
	v[665] = v[161] * v[625] + v[160] * v[629];
	v[666] = v[161] * v[626] + v[160] * v[630];
	v[667] = v[158] * v[623] + v[157] * v[627];
	v[1034] = v[281] * v[667];
	v[668] = v[158] * v[624] + v[157] * v[628];
	v[1086] = v[281] * v[668];
	v[669] = v[158] * v[625] + v[157] * v[629];
	v[1138] = v[281] * v[669];
	v[670] = v[158] * v[626] + v[157] * v[630];
	v[1190] = v[281] * v[670];
	v[671] = v[164] * v[607] + v[163] * v[611];
	v[1030] = v[281] * v[671];
	v[672] = v[164] * v[608] + v[163] * v[612];
	v[1082] = v[281] * v[672];
	v[673] = v[164] * v[609] + v[163] * v[613];
	v[1134] = v[281] * v[673];
	v[674] = v[164] * v[610] + v[163] * v[614];
	v[1186] = v[281] * v[674];
	v[675] = v[161] * v[607] + v[160] * v[611];
	v[1035] = v[281] * v[675];
	v[676] = v[161] * v[608] + v[160] * v[612];
	v[1087] = v[281] * v[676];
	v[677] = v[161] * v[609] + v[160] * v[613];
	v[1139] = v[281] * v[677];
	v[678] = v[161] * v[610] + v[160] * v[614];
	v[1191] = v[281] * v[678];
	v[679] = v[158] * v[607] + v[157] * v[611];
	v[680] = v[158] * v[608] + v[157] * v[612];
	v[681] = v[158] * v[609] + v[157] * v[613];
	v[682] = v[158] * v[610] + v[157] * v[614];
	v[683] = v[155] * v[631] + v[154] * v[635];
	v[684] = v[155] * v[632] + v[154] * v[636];
	v[685] = v[155] * v[633] + v[154] * v[637];
	v[686] = v[155] * v[634] + v[154] * v[638];
	v[687] = v[152] * v[631] + v[151] * v[635];
	v[1013] = v[262] * v[687];
	v[688] = v[152] * v[632] + v[151] * v[636];
	v[1065] = v[262] * v[688];
	v[689] = v[152] * v[633] + v[151] * v[637];
	v[1117] = v[262] * v[689];
	v[690] = v[152] * v[634] + v[151] * v[638];
	v[1169] = v[262] * v[690];
	v[691] = v[149] * v[631] + v[148] * v[635];
	v[1016] = v[262] * v[691];
	v[692] = v[149] * v[632] + v[148] * v[636];
	v[1068] = v[262] * v[692];
	v[693] = v[149] * v[633] + v[148] * v[637];
	v[1120] = v[262] * v[693];
	v[694] = v[149] * v[634] + v[148] * v[638];
	v[1172] = v[262] * v[694];
	v[695] = v[155] * v[615] + v[154] * v[619];
	v[1014] = v[262] * v[695];
	v[696] = v[155] * v[616] + v[154] * v[620];
	v[1066] = v[262] * v[696];
	v[697] = v[155] * v[617] + v[154] * v[621];
	v[1118] = v[262] * v[697];
	v[698] = v[155] * v[618] + v[154] * v[622];
	v[1170] = v[262] * v[698];
	v[699] = v[152] * v[615] + v[151] * v[619];
	v[700] = v[152] * v[616] + v[151] * v[620];
	v[701] = v[152] * v[617] + v[151] * v[621];
	v[702] = v[152] * v[618] + v[151] * v[622];
	v[703] = v[149] * v[615] + v[148] * v[619];
	v[1021] = v[262] * v[703];
	v[704] = v[149] * v[616] + v[148] * v[620];
	v[1073] = v[262] * v[704];
	v[705] = v[149] * v[617] + v[148] * v[621];
	v[1125] = v[262] * v[705];
	v[706] = v[149] * v[618] + v[148] * v[622];
	v[1177] = v[262] * v[706];
	v[707] = v[155] * v[599] + v[154] * v[603];
	v[1017] = v[262] * v[707];
	v[708] = v[155] * v[600] + v[154] * v[604];
	v[1069] = v[262] * v[708];
	v[709] = v[155] * v[601] + v[154] * v[605];
	v[1121] = v[262] * v[709];
	v[710] = v[155] * v[602] + v[154] * v[606];
	v[1173] = v[262] * v[710];
	v[711] = v[152] * v[599] + v[151] * v[603];
	v[1022] = v[262] * v[711];
	v[712] = v[152] * v[600] + v[151] * v[604];
	v[1074] = v[262] * v[712];
	v[713] = v[152] * v[601] + v[151] * v[605];
	v[1126] = v[262] * v[713];
	v[714] = v[152] * v[602] + v[151] * v[606];
	v[1178] = v[262] * v[714];
	v[715] = v[149] * v[599] + v[148] * v[603];
	v[716] = v[149] * v[600] + v[148] * v[604];
	v[717] = v[149] * v[601] + v[148] * v[605];
	v[718] = v[149] * v[602] + v[148] * v[606];
	v[1193] = -(v[4741] * v[682]);
	v[1192] = -(v[4741] * v[666]);
	v[1141] = -(v[4741] * v[681]);
	v[1140] = -(v[4741] * v[665]);
	v[1089] = -(v[4741] * v[680]);
	v[1088] = -(v[4741] * v[664]);
	v[1037] = -(v[4741] * v[679]);
	v[1036] = -(v[4741] * v[663]);
	v[725] = v[1026] + v[1027];
	v[726] = v[1078] + v[1079];
	v[727] = v[1130] + v[1131];
	v[728] = v[1182] + v[1183];
	v[733] = v[1029] + v[1030];
	v[734] = v[1081] + v[1082];
	v[735] = v[1133] + v[1134];
	v[736] = v[1185] + v[1186];
	v[738] = v[1034] + v[1035];
	v[739] = v[1086] + v[1087];
	v[740] = v[1138] + v[1139];
	v[741] = v[1190] + v[1191];
	v[1188] = (v[4784] * v[650] + v[4785] * v[666] + v[4786] * v[682] + v[654] * v[722] + v[658] * v[723] + v[662] * v[724]
		+ v[670] * v[731] + v[674] * v[732] + v[678] * v[737]) * v[749];
	v[4817] = v[1188] - v[4741] * v[650];
	v[1136] = (v[4784] * v[649] + v[4785] * v[665] + v[4786] * v[681] + v[653] * v[722] + v[657] * v[723] + v[661] * v[724]
		+ v[669] * v[731] + v[673] * v[732] + v[677] * v[737]) * v[749];
	v[4813] = v[1136] - v[4741] * v[649];
	v[1084] = (v[4784] * v[648] + v[4785] * v[664] + v[4786] * v[680] + v[652] * v[722] + v[656] * v[723] + v[660] * v[724]
		+ v[668] * v[731] + v[672] * v[732] + v[676] * v[737]) * v[749];
	v[4809] = v[1084] - v[4741] * v[648];
	v[1032] = (v[4784] * v[647] + v[4785] * v[663] + v[4786] * v[679] + v[651] * v[722] + v[655] * v[723] + v[659] * v[724]
		+ v[667] * v[731] + v[671] * v[732] + v[675] * v[737]) * v[749];
	v[4805] = v[1032] - v[4741] * v[647];
	v[1180] = -(v[4737] * v[718]);
	v[1179] = -(v[4737] * v[702]);
	v[1128] = -(v[4737] * v[717]);
	v[1127] = -(v[4737] * v[701]);
	v[1076] = -(v[4737] * v[716]);
	v[1075] = -(v[4737] * v[700]);
	v[1024] = -(v[4737] * v[715]);
	v[1023] = -(v[4737] * v[699]);
	v[756] = v[1013] + v[1014];
	v[757] = v[1065] + v[1066];
	v[758] = v[1117] + v[1118];
	v[759] = v[1169] + v[1170];
	v[764] = v[1016] + v[1017];
	v[765] = v[1068] + v[1069];
	v[766] = v[1120] + v[1121];
	v[767] = v[1172] + v[1173];
	v[769] = v[1021] + v[1022];
	v[770] = v[1073] + v[1074];
	v[771] = v[1125] + v[1126];
	v[772] = v[1177] + v[1178];
	v[1175] = (v[4787] * v[686] + v[4788] * v[702] + v[4789] * v[718] + v[690] * v[753] + v[694] * v[754] + v[698] * v[755]
		+ v[706] * v[762] + v[710] * v[763] + v[714] * v[768]) * v[780];
	v[4816] = v[1175] - v[4737] * v[686];
	v[1123] = (v[4787] * v[685] + v[4788] * v[701] + v[4789] * v[717] + v[689] * v[753] + v[693] * v[754] + v[697] * v[755]
		+ v[705] * v[762] + v[709] * v[763] + v[713] * v[768]) * v[780];
	v[4812] = v[1123] - v[4737] * v[685];
	v[1071] = (v[4787] * v[684] + v[4788] * v[700] + v[4789] * v[716] + v[688] * v[753] + v[692] * v[754] + v[696] * v[755]
		+ v[704] * v[762] + v[708] * v[763] + v[712] * v[768]) * v[780];
	v[4808] = v[1071] - v[4737] * v[684];
	v[1019] = (v[4787] * v[683] + v[4788] * v[699] + v[4789] * v[715] + v[687] * v[753] + v[691] * v[754] + v[695] * v[755]
		+ v[703] * v[762] + v[707] * v[763] + v[711] * v[768]) * v[780];
	v[4804] = v[1019] - v[4737] * v[683];
	v[781] = v[252] * v[591];
	v[782] = v[4790] * v[509] + v[252] * v[592];
	v[783] = v[252] * v[593];
	v[784] = v[252] * v[594];
	v[785] = v[250] * v[591];
	v[786] = v[4790] * v[508] + v[250] * v[592];
	v[787] = v[250] * v[593];
	v[788] = v[250] * v[594];
	v[789] = v[252] * v[595];
	v[790] = v[4791] * v[509] + v[252] * v[596];
	v[791] = v[252] * v[597];
	v[792] = v[252] * v[598];
	v[793] = v[250] * v[595];
	v[794] = v[4791] * v[508] + v[250] * v[596];
	v[795] = v[250] * v[597];
	v[796] = v[250] * v[598];
	v[797] = v[252] * v[574];
	v[798] = v[4792] * v[509] + v[252] * v[575];
	v[799] = v[252] * v[576];
	v[800] = v[252] * v[577];
	v[801] = v[250] * v[574];
	v[802] = v[4792] * v[508] + v[250] * v[575];
	v[803] = v[250] * v[576];
	v[804] = v[250] * v[577];
	v[805] = v[252] * v[578];
	v[806] = v[4793] * v[509] + v[252] * v[579];
	v[807] = v[252] * v[580];
	v[808] = v[252] * v[581];
	v[809] = v[250] * v[578];
	v[810] = v[4793] * v[508] + v[250] * v[579];
	v[811] = v[250] * v[580];
	v[812] = v[250] * v[581];
	v[813] = v[252] * v[557];
	v[814] = v[4794] * v[509] + v[252] * v[558];
	v[815] = v[252] * v[559];
	v[816] = v[252] * v[560];
	v[817] = v[250] * v[557];
	v[818] = v[4794] * v[508] + v[250] * v[558];
	v[819] = v[250] * v[559];
	v[820] = v[250] * v[560];
	v[821] = v[252] * v[561];
	v[822] = v[4795] * v[509] + v[252] * v[562];
	v[823] = v[252] * v[563];
	v[824] = v[252] * v[564];
	v[825] = v[250] * v[561];
	v[826] = v[4795] * v[508] + v[250] * v[562];
	v[827] = v[250] * v[563];
	v[828] = v[250] * v[564];
	v[829] = v[101] * v[821] + v[100] * v[825];
	v[830] = v[101] * v[822] + v[100] * v[826];
	v[831] = v[101] * v[823] + v[100] * v[827];
	v[832] = v[101] * v[824] + v[100] * v[828];
	v[833] = v[825] * v[97] + v[821] * v[98];
	v[1000] = v[194] * v[833];
	v[834] = v[826] * v[97] + v[822] * v[98];
	v[1052] = v[194] * v[834];
	v[835] = v[827] * v[97] + v[823] * v[98];
	v[1104] = v[194] * v[835];
	v[836] = v[828] * v[97] + v[824] * v[98];
	v[1156] = v[194] * v[836];
	v[837] = v[825] * v[94] + v[821] * v[95];
	v[1003] = v[194] * v[837];
	v[838] = v[826] * v[94] + v[822] * v[95];
	v[1055] = v[194] * v[838];
	v[839] = v[827] * v[94] + v[823] * v[95];
	v[1107] = v[194] * v[839];
	v[840] = v[828] * v[94] + v[824] * v[95];
	v[1159] = v[194] * v[840];
	v[841] = v[101] * v[805] + v[100] * v[809];
	v[1001] = v[194] * v[841];
	v[842] = v[101] * v[806] + v[100] * v[810];
	v[1053] = v[194] * v[842];
	v[843] = v[101] * v[807] + v[100] * v[811];
	v[1105] = v[194] * v[843];
	v[844] = v[101] * v[808] + v[100] * v[812];
	v[1157] = v[194] * v[844];
	v[845] = v[809] * v[97] + v[805] * v[98];
	v[846] = v[810] * v[97] + v[806] * v[98];
	v[847] = v[811] * v[97] + v[807] * v[98];
	v[848] = v[812] * v[97] + v[808] * v[98];
	v[849] = v[809] * v[94] + v[805] * v[95];
	v[1008] = v[194] * v[849];
	v[850] = v[810] * v[94] + v[806] * v[95];
	v[1060] = v[194] * v[850];
	v[851] = v[811] * v[94] + v[807] * v[95];
	v[1112] = v[194] * v[851];
	v[852] = v[812] * v[94] + v[808] * v[95];
	v[1164] = v[194] * v[852];
	v[853] = v[101] * v[789] + v[100] * v[793];
	v[1004] = v[194] * v[853];
	v[854] = v[101] * v[790] + v[100] * v[794];
	v[1056] = v[194] * v[854];
	v[855] = v[101] * v[791] + v[100] * v[795];
	v[1108] = v[194] * v[855];
	v[856] = v[101] * v[792] + v[100] * v[796];
	v[1160] = v[194] * v[856];
	v[857] = v[793] * v[97] + v[789] * v[98];
	v[1009] = v[194] * v[857];
	v[858] = v[794] * v[97] + v[790] * v[98];
	v[1061] = v[194] * v[858];
	v[859] = v[795] * v[97] + v[791] * v[98];
	v[1113] = v[194] * v[859];
	v[860] = v[796] * v[97] + v[792] * v[98];
	v[1165] = v[194] * v[860];
	v[861] = v[793] * v[94] + v[789] * v[95];
	v[862] = v[794] * v[94] + v[790] * v[95];
	v[863] = v[795] * v[94] + v[791] * v[95];
	v[864] = v[796] * v[94] + v[792] * v[95];
	v[865] = v[817] * v[91] + v[813] * v[92];
	v[866] = v[818] * v[91] + v[814] * v[92];
	v[867] = v[819] * v[91] + v[815] * v[92];
	v[868] = v[820] * v[91] + v[816] * v[92];
	v[869] = v[817] * v[88] + v[813] * v[89];
	v[987] = v[175] * v[869];
	v[870] = v[818] * v[88] + v[814] * v[89];
	v[1039] = v[175] * v[870];
	v[871] = v[819] * v[88] + v[815] * v[89];
	v[1091] = v[175] * v[871];
	v[872] = v[820] * v[88] + v[816] * v[89];
	v[1143] = v[175] * v[872];
	v[873] = v[817] * v[85] + v[813] * v[86];
	v[990] = v[175] * v[873];
	v[874] = v[818] * v[85] + v[814] * v[86];
	v[1042] = v[175] * v[874];
	v[875] = v[819] * v[85] + v[815] * v[86];
	v[1094] = v[175] * v[875];
	v[876] = v[820] * v[85] + v[816] * v[86];
	v[1146] = v[175] * v[876];
	v[877] = v[801] * v[91] + v[797] * v[92];
	v[988] = v[175] * v[877];
	v[878] = v[802] * v[91] + v[798] * v[92];
	v[1040] = v[175] * v[878];
	v[879] = v[803] * v[91] + v[799] * v[92];
	v[1092] = v[175] * v[879];
	v[880] = v[804] * v[91] + v[800] * v[92];
	v[1144] = v[175] * v[880];
	v[881] = v[801] * v[88] + v[797] * v[89];
	v[882] = v[802] * v[88] + v[798] * v[89];
	v[883] = v[803] * v[88] + v[799] * v[89];
	v[884] = v[804] * v[88] + v[800] * v[89];
	v[885] = v[801] * v[85] + v[797] * v[86];
	v[995] = v[175] * v[885];
	v[886] = v[802] * v[85] + v[798] * v[86];
	v[1047] = v[175] * v[886];
	v[887] = v[803] * v[85] + v[799] * v[86];
	v[1099] = v[175] * v[887];
	v[888] = v[804] * v[85] + v[800] * v[86];
	v[1151] = v[175] * v[888];
	v[889] = v[785] * v[91] + v[781] * v[92];
	v[991] = v[175] * v[889];
	v[890] = v[786] * v[91] + v[782] * v[92];
	v[1043] = v[175] * v[890];
	v[891] = v[787] * v[91] + v[783] * v[92];
	v[1095] = v[175] * v[891];
	v[892] = v[788] * v[91] + v[784] * v[92];
	v[1147] = v[175] * v[892];
	v[893] = v[785] * v[88] + v[781] * v[89];
	v[996] = v[175] * v[893];
	v[894] = v[786] * v[88] + v[782] * v[89];
	v[1048] = v[175] * v[894];
	v[895] = v[787] * v[88] + v[783] * v[89];
	v[1100] = v[175] * v[895];
	v[896] = v[788] * v[88] + v[784] * v[89];
	v[1152] = v[175] * v[896];
	v[897] = v[785] * v[85] + v[781] * v[86];
	v[898] = v[786] * v[85] + v[782] * v[86];
	v[899] = v[787] * v[85] + v[783] * v[86];
	v[900] = v[788] * v[85] + v[784] * v[86];
	v[1167] = -(v[4728] * v[864]);
	v[1166] = -(v[4728] * v[848]);
	v[1115] = -(v[4728] * v[863]);
	v[1114] = -(v[4728] * v[847]);
	v[1063] = -(v[4728] * v[862]);
	v[1062] = -(v[4728] * v[846]);
	v[1011] = -(v[4728] * v[861]);
	v[1010] = -(v[4728] * v[845]);
	v[907] = v[1000] + v[1001];
	v[908] = v[1052] + v[1053];
	v[909] = v[1104] + v[1105];
	v[910] = v[1156] + v[1157];
	v[915] = v[1003] + v[1004];
	v[916] = v[1055] + v[1056];
	v[917] = v[1107] + v[1108];
	v[918] = v[1159] + v[1160];
	v[920] = v[1008] + v[1009];
	v[921] = v[1060] + v[1061];
	v[922] = v[1112] + v[1113];
	v[923] = v[1164] + v[1165];
	v[1162] = (v[4796] * v[832] + v[4797] * v[848] + v[4798] * v[864] + v[836] * v[904] + v[840] * v[905] + v[844] * v[906]
		+ v[852] * v[913] + v[856] * v[914] + v[860] * v[919]) * v[931];
	v[4815] = v[1162] - v[4728] * v[832];
	v[1110] = (v[4796] * v[831] + v[4797] * v[847] + v[4798] * v[863] + v[835] * v[904] + v[839] * v[905] + v[843] * v[906]
		+ v[851] * v[913] + v[855] * v[914] + v[859] * v[919]) * v[931];
	v[4811] = v[1110] - v[4728] * v[831];
	v[1058] = (v[4796] * v[830] + v[4797] * v[846] + v[4798] * v[862] + v[834] * v[904] + v[838] * v[905] + v[842] * v[906]
		+ v[850] * v[913] + v[854] * v[914] + v[858] * v[919]) * v[931];
	v[4807] = v[1058] - v[4728] * v[830];
	v[1006] = (v[4796] * v[829] + v[4797] * v[845] + v[4798] * v[861] + v[833] * v[904] + v[837] * v[905] + v[841] * v[906]
		+ v[849] * v[913] + v[853] * v[914] + v[857] * v[919]) * v[931];
	v[4803] = v[1006] - v[4728] * v[829];
	v[1154] = -(v[4724] * v[900]);
	v[1153] = -(v[4724] * v[884]);
	v[1102] = -(v[4724] * v[899]);
	v[1101] = -(v[4724] * v[883]);
	v[1050] = -(v[4724] * v[898]);
	v[1049] = -(v[4724] * v[882]);
	v[998] = -(v[4724] * v[897]);
	v[997] = -(v[4724] * v[881]);
	v[938] = v[987] + v[988];
	v[939] = v[1039] + v[1040];
	v[940] = v[1091] + v[1092];
	v[941] = v[1143] + v[1144];
	v[946] = v[990] + v[991];
	v[947] = v[1042] + v[1043];
	v[948] = v[1094] + v[1095];
	v[949] = v[1146] + v[1147];
	v[951] = v[995] + v[996];
	v[952] = v[1047] + v[1048];
	v[953] = v[1099] + v[1100];
	v[954] = v[1151] + v[1152];
	v[1149] = (v[4799] * v[868] + v[4800] * v[884] + v[4801] * v[900] + v[872] * v[935] + v[876] * v[936] + v[880] * v[937]
		+ v[888] * v[944] + v[892] * v[945] + v[896] * v[950]) * v[962];
	v[4814] = v[1149] - v[4724] * v[868];
	v[1097] = (v[4799] * v[867] + v[4800] * v[883] + v[4801] * v[899] + v[871] * v[935] + v[875] * v[936] + v[879] * v[937]
		+ v[887] * v[944] + v[891] * v[945] + v[895] * v[950]) * v[962];
	v[4810] = v[1097] - v[4724] * v[867];
	v[1045] = (v[4799] * v[866] + v[4800] * v[882] + v[4801] * v[898] + v[870] * v[935] + v[874] * v[936] + v[878] * v[937]
		+ v[886] * v[944] + v[890] * v[945] + v[894] * v[950]) * v[962];
	v[4806] = v[1045] - v[4724] * v[866];
	v[993] = (v[4799] * v[865] + v[4800] * v[881] + v[4801] * v[897] + v[869] * v[935] + v[873] * v[936] + v[877] * v[937]
		+ v[885] * v[944] + v[889] * v[945] + v[893] * v[950]) * v[962];
	v[4802] = -(v[4724] * v[865]) + v[993];
	v[989] = v[946] * v[983] + v[951] * v[984] + v[987] - v[988] + v[986] * (v[4802] + v[997]);
	v[994] = v[951] * v[982] + v[938] * v[983] - v[990] + v[991] + v[985] * (v[4802] + v[998]);
	v[999] = v[946] * v[982] + v[938] * v[984] + v[995] - v[996] + v[981] * (v[993] + v[997] + v[998]);
	v[1002] = v[1000] - v[1001] + v[915] * v[977] + v[920] * v[978] + (v[1010] + v[4803]) * v[980];
	v[1007] = -v[1003] + v[1004] + v[920] * v[976] + v[907] * v[977] + (v[1011] + v[4803]) * v[979];
	v[1012] = v[1008] - v[1009] + (v[1006] + v[1010] + v[1011]) * v[975] + v[915] * v[976] + v[907] * v[978];
	v[1015] = v[1013] - v[1014] + v[764] * v[971] + v[769] * v[972] + (v[1023] + v[4804]) * v[974];
	v[1020] = -v[1016] + v[1017] + v[769] * v[970] + v[756] * v[971] + (v[1024] + v[4804]) * v[973];
	v[1025] = v[1021] - v[1022] + (v[1019] + v[1023] + v[1024]) * v[969] + v[764] * v[970] + v[756] * v[972];
	v[1028] = v[1026] - v[1027] + v[733] * v[965] + v[738] * v[966] + (v[1036] + v[4805]) * v[968];
	v[1033] = -v[1029] + v[1030] + v[738] * v[964] + v[725] * v[965] + (v[1037] + v[4805]) * v[967];
	v[1038] = v[1034] - v[1035] + (v[1032] + v[1036] + v[1037]) * v[963] + v[733] * v[964] + v[725] * v[966];
	v[1041] = v[1039] - v[1040] + v[947] * v[983] + v[952] * v[984] + (v[1049] + v[4806]) * v[986];
	v[1046] = -v[1042] + v[1043] + v[952] * v[982] + v[939] * v[983] + (v[1050] + v[4806]) * v[985];
	v[1051] = v[1047] - v[1048] + (v[1045] + v[1049] + v[1050]) * v[981] + v[947] * v[982] + v[939] * v[984];
	v[1054] = v[1052] - v[1053] + v[916] * v[977] + v[921] * v[978] + (v[1062] + v[4807]) * v[980];
	v[1059] = -v[1055] + v[1056] + v[921] * v[976] + v[908] * v[977] + (v[1063] + v[4807]) * v[979];
	v[1064] = v[1060] - v[1061] + (v[1058] + v[1062] + v[1063]) * v[975] + v[916] * v[976] + v[908] * v[978];
	v[1067] = v[1065] - v[1066] + v[765] * v[971] + v[770] * v[972] + (v[1075] + v[4808]) * v[974];
	v[1072] = -v[1068] + v[1069] + v[770] * v[970] + v[757] * v[971] + (v[1076] + v[4808]) * v[973];
	v[1077] = v[1073] - v[1074] + (v[1071] + v[1075] + v[1076]) * v[969] + v[765] * v[970] + v[757] * v[972];
	v[1080] = v[1078] - v[1079] + v[734] * v[965] + v[739] * v[966] + (v[1088] + v[4809]) * v[968];
	v[1085] = -v[1081] + v[1082] + v[739] * v[964] + v[726] * v[965] + (v[1089] + v[4809]) * v[967];
	v[1090] = v[1086] - v[1087] + (v[1084] + v[1088] + v[1089]) * v[963] + v[734] * v[964] + v[726] * v[966];
	v[1093] = v[1091] - v[1092] + v[948] * v[983] + v[953] * v[984] + (v[1101] + v[4810]) * v[986];
	v[1098] = -v[1094] + v[1095] + v[953] * v[982] + v[940] * v[983] + (v[1102] + v[4810]) * v[985];
	v[1103] = v[1099] - v[1100] + (v[1097] + v[1101] + v[1102]) * v[981] + v[948] * v[982] + v[940] * v[984];
	v[1106] = v[1104] - v[1105] + v[917] * v[977] + v[922] * v[978] + (v[1114] + v[4811]) * v[980];
	v[1111] = -v[1107] + v[1108] + v[922] * v[976] + v[909] * v[977] + (v[1115] + v[4811]) * v[979];
	v[1116] = v[1112] - v[1113] + (v[1110] + v[1114] + v[1115]) * v[975] + v[917] * v[976] + v[909] * v[978];
	v[1119] = v[1117] - v[1118] + v[766] * v[971] + v[771] * v[972] + (v[1127] + v[4812]) * v[974];
	v[1124] = -v[1120] + v[1121] + v[771] * v[970] + v[758] * v[971] + (v[1128] + v[4812]) * v[973];
	v[1129] = v[1125] - v[1126] + (v[1123] + v[1127] + v[1128]) * v[969] + v[766] * v[970] + v[758] * v[972];
	v[1132] = v[1130] - v[1131] + v[735] * v[965] + v[740] * v[966] + (v[1140] + v[4813]) * v[968];
	v[1137] = -v[1133] + v[1134] + v[740] * v[964] + v[727] * v[965] + (v[1141] + v[4813]) * v[967];
	v[1142] = v[1138] - v[1139] + (v[1136] + v[1140] + v[1141]) * v[963] + v[735] * v[964] + v[727] * v[966];
	v[1145] = v[1143] - v[1144] + v[949] * v[983] + v[954] * v[984] + (v[1153] + v[4814]) * v[986];
	v[1150] = -v[1146] + v[1147] + v[954] * v[982] + v[941] * v[983] + (v[1154] + v[4814]) * v[985];
	v[1155] = v[1151] - v[1152] + (v[1149] + v[1153] + v[1154]) * v[981] + v[949] * v[982] + v[941] * v[984];
	v[1158] = v[1156] - v[1157] + v[918] * v[977] + v[923] * v[978] + (v[1166] + v[4815]) * v[980];
	v[1163] = -v[1159] + v[1160] + v[923] * v[976] + v[910] * v[977] + (v[1167] + v[4815]) * v[979];
	v[1168] = v[1164] - v[1165] + (v[1162] + v[1166] + v[1167]) * v[975] + v[918] * v[976] + v[910] * v[978];
	v[1171] = v[1169] - v[1170] + v[767] * v[971] + v[772] * v[972] + (v[1179] + v[4816]) * v[974];
	v[1176] = -v[1172] + v[1173] + v[772] * v[970] + v[759] * v[971] + (v[1180] + v[4816]) * v[973];
	v[1181] = v[1177] - v[1178] + (v[1175] + v[1179] + v[1180]) * v[969] + v[767] * v[970] + v[759] * v[972];
	v[1184] = v[1182] - v[1183] + v[736] * v[965] + v[741] * v[966] + (v[1192] + v[4817]) * v[968];
	v[1189] = -v[1185] + v[1186] + v[741] * v[964] + v[728] * v[965] + (v[1193] + v[4817]) * v[967];
	v[1194] = v[1190] - v[1191] + (v[1188] + v[1192] + v[1193]) * v[963] + v[736] * v[964] + v[728] * v[966];
	v[1211] = -(v[22] * v[591]) - v[23] * v[592] - v[24] * v[593] - v[25] * v[594];
	v[1212] = -(v[22] * v[574]) - v[23] * v[575] - v[24] * v[576] - v[25] * v[577];
	v[1213] = -(v[22] * v[557]) - v[23] * v[558] - v[24] * v[559] - v[25] * v[560];
	v[1214] = -(v[1041] * v[23]) - v[1093] * v[24] - v[1145] * v[25] - v[22] * v[989];
	v[1215] = -(v[1046] * v[23]) - v[1098] * v[24] - v[1150] * v[25] - v[22] * v[994];
	v[1216] = -(v[1051] * v[23]) - v[1103] * v[24] - v[1155] * v[25] - v[22] * v[999];
	v[1217] = -(v[22] * v[595]) - v[23] * v[596] - v[24] * v[597] - v[25] * v[598];
	v[1218] = -(v[22] * v[578]) - v[23] * v[579] - v[24] * v[580] - v[25] * v[581];
	v[1219] = -(v[22] * v[561]) - v[23] * v[562] - v[24] * v[563] - v[25] * v[564];
	v[1220] = -(v[1002] * v[22]) - v[1054] * v[23] - v[1106] * v[24] - v[1158] * v[25];
	v[1221] = -(v[1007] * v[22]) - v[1059] * v[23] - v[1111] * v[24] - v[1163] * v[25];
	v[1222] = -(v[1012] * v[22]) - v[1064] * v[23] - v[1116] * v[24] - v[1168] * v[25];
	v[1223] = -(v[22] * v[582]) - v[23] * v[583] - v[24] * v[584] - v[25] * v[585];
	v[1224] = -(v[22] * v[565]) - v[23] * v[566] - v[24] * v[567] - v[25] * v[568];
	v[1225] = -(v[22] * v[548]) - v[23] * v[549] - v[24] * v[550] - v[25] * v[551];
	v[1226] = -(v[1015] * v[22]) - v[1067] * v[23] - v[1119] * v[24] - v[1171] * v[25];
	v[1227] = -(v[1020] * v[22]) - v[1072] * v[23] - v[1124] * v[24] - v[1176] * v[25];
	v[1228] = -(v[1025] * v[22]) - v[1077] * v[23] - v[1129] * v[24] - v[1181] * v[25];
	v[1229] = -(v[22] * v[586]) - v[23] * v[587] - v[24] * v[589] - v[25] * v[590];
	v[1230] = -(v[22] * v[569]) - v[23] * v[570] - v[24] * v[572] - v[25] * v[573];
	v[1231] = -(v[22] * v[552]) - v[23] * v[553] - v[24] * v[555] - v[25] * v[556];
	v[1232] = -(v[1028] * v[22]) - v[1080] * v[23] - v[1132] * v[24] - v[1184] * v[25];
	v[1233] = -(v[1033] * v[22]) - v[1085] * v[23] - v[1137] * v[24] - v[1189] * v[25];
	v[1234] = -(v[1038] * v[22]) - v[1090] * v[23] - v[1142] * v[24] - v[1194] * v[25];
	v[9878] = v[1211];
	v[9879] = v[1212];
	v[9880] = v[1213];
	v[9881] = v[1214];
	v[9882] = v[1215];
	v[9883] = v[1216];
	v[9884] = v[1217];
	v[9885] = v[1218];
	v[9886] = v[1219];
	v[9887] = v[1220];
	v[9888] = v[1221];
	v[9889] = v[1222];
	v[9890] = v[1223];
	v[9891] = v[1224];
	v[9892] = v[1225];
	v[9893] = v[1226];
	v[9894] = v[1227];
	v[9895] = v[1228];
	v[9896] = v[1229];
	v[9897] = v[1230];
	v[9898] = v[1231];
	v[9899] = v[1232];
	v[9900] = v[1233];
	v[9901] = v[1234];
	v[1732] = -(v[1211] * v[349]) - v[1212] * v[350] - v[1213] * v[351] - v[1217] * v[352] - v[1218] * v[353]
		- v[1219] * v[354] - v[1223] * v[355] - v[1224] * v[356] - v[1225] * v[357] - v[1229] * v[358] - v[1230] * v[359]
		- v[1231] * v[360] - v[1214] * v[366] - v[1215] * v[376] - v[1216] * v[386] - v[1220] * v[392] - v[1221] * v[402]
		- v[1222] * v[412] - v[1226] * v[418] - v[1227] * v[428] - v[1228] * v[438] - v[1232] * v[444] - v[1233] * v[454]
		- v[1234] * v[464];
	v[1235] = -(v[26] * v[591]) - v[27] * v[592] - v[28] * v[593] - v[29] * v[594];
	v[1236] = -(v[26] * v[574]) - v[27] * v[575] - v[28] * v[576] - v[29] * v[577];
	v[1237] = -(v[26] * v[557]) - v[27] * v[558] - v[28] * v[559] - v[29] * v[560];
	v[1238] = -(v[1041] * v[27]) - v[1093] * v[28] - v[1145] * v[29] - v[26] * v[989];
	v[1239] = -(v[1046] * v[27]) - v[1098] * v[28] - v[1150] * v[29] - v[26] * v[994];
	v[1240] = -(v[1051] * v[27]) - v[1103] * v[28] - v[1155] * v[29] - v[26] * v[999];
	v[1241] = -(v[26] * v[595]) - v[27] * v[596] - v[28] * v[597] - v[29] * v[598];
	v[1242] = -(v[26] * v[578]) - v[27] * v[579] - v[28] * v[580] - v[29] * v[581];
	v[1243] = -(v[26] * v[561]) - v[27] * v[562] - v[28] * v[563] - v[29] * v[564];
	v[1244] = -(v[1002] * v[26]) - v[1054] * v[27] - v[1106] * v[28] - v[1158] * v[29];
	v[1245] = -(v[1007] * v[26]) - v[1059] * v[27] - v[1111] * v[28] - v[1163] * v[29];
	v[1246] = -(v[1012] * v[26]) - v[1064] * v[27] - v[1116] * v[28] - v[1168] * v[29];
	v[1247] = -(v[26] * v[582]) - v[27] * v[583] - v[28] * v[584] - v[29] * v[585];
	v[1248] = -(v[26] * v[565]) - v[27] * v[566] - v[28] * v[567] - v[29] * v[568];
	v[1249] = -(v[26] * v[548]) - v[27] * v[549] - v[28] * v[550] - v[29] * v[551];
	v[1250] = -(v[1015] * v[26]) - v[1067] * v[27] - v[1119] * v[28] - v[1171] * v[29];
	v[1251] = -(v[1020] * v[26]) - v[1072] * v[27] - v[1124] * v[28] - v[1176] * v[29];
	v[1252] = -(v[1025] * v[26]) - v[1077] * v[27] - v[1129] * v[28] - v[1181] * v[29];
	v[1253] = -(v[26] * v[586]) - v[27] * v[587] - v[28] * v[589] - v[29] * v[590];
	v[1254] = -(v[26] * v[569]) - v[27] * v[570] - v[28] * v[572] - v[29] * v[573];
	v[1255] = -(v[26] * v[552]) - v[27] * v[553] - v[28] * v[555] - v[29] * v[556];
	v[1256] = -(v[1028] * v[26]) - v[1080] * v[27] - v[1132] * v[28] - v[1184] * v[29];
	v[1257] = -(v[1033] * v[26]) - v[1085] * v[27] - v[1137] * v[28] - v[1189] * v[29];
	v[1258] = -(v[1038] * v[26]) - v[1090] * v[27] - v[1142] * v[28] - v[1194] * v[29];
	v[9830] = v[1235];
	v[9831] = v[1236];
	v[9832] = v[1237];
	v[9833] = v[1238];
	v[9834] = v[1239];
	v[9835] = v[1240];
	v[9836] = v[1241];
	v[9837] = v[1242];
	v[9838] = v[1243];
	v[9839] = v[1244];
	v[9840] = v[1245];
	v[9841] = v[1246];
	v[9842] = v[1247];
	v[9843] = v[1248];
	v[9844] = v[1249];
	v[9845] = v[1250];
	v[9846] = v[1251];
	v[9847] = v[1252];
	v[9848] = v[1253];
	v[9849] = v[1254];
	v[9850] = v[1255];
	v[9851] = v[1256];
	v[9852] = v[1257];
	v[9853] = v[1258];
	v[1734] = -(v[1235] * v[349]) - v[1236] * v[350] - v[1237] * v[351] - v[1241] * v[352] - v[1242] * v[353]
		- v[1243] * v[354] - v[1247] * v[355] - v[1248] * v[356] - v[1249] * v[357] - v[1253] * v[358] - v[1254] * v[359]
		- v[1255] * v[360] - v[1238] * v[366] - v[1239] * v[376] - v[1240] * v[386] - v[1244] * v[392] - v[1245] * v[402]
		- v[1246] * v[412] - v[1250] * v[418] - v[1251] * v[428] - v[1252] * v[438] - v[1256] * v[444] - v[1257] * v[454]
		- v[1258] * v[464];
	v[1259] = -(v[30] * v[591]) - v[31] * v[592] - v[32] * v[593] - v[33] * v[594];
	v[1260] = -(v[30] * v[574]) - v[31] * v[575] - v[32] * v[576] - v[33] * v[577];
	v[1261] = -(v[30] * v[557]) - v[31] * v[558] - v[32] * v[559] - v[33] * v[560];
	v[1262] = -(v[1041] * v[31]) - v[1093] * v[32] - v[1145] * v[33] - v[30] * v[989];
	v[1263] = -(v[1046] * v[31]) - v[1098] * v[32] - v[1150] * v[33] - v[30] * v[994];
	v[1264] = -(v[1051] * v[31]) - v[1103] * v[32] - v[1155] * v[33] - v[30] * v[999];
	v[1265] = -(v[30] * v[595]) - v[31] * v[596] - v[32] * v[597] - v[33] * v[598];
	v[1266] = -(v[30] * v[578]) - v[31] * v[579] - v[32] * v[580] - v[33] * v[581];
	v[1267] = -(v[30] * v[561]) - v[31] * v[562] - v[32] * v[563] - v[33] * v[564];
	v[1268] = -(v[1002] * v[30]) - v[1054] * v[31] - v[1106] * v[32] - v[1158] * v[33];
	v[1269] = -(v[1007] * v[30]) - v[1059] * v[31] - v[1111] * v[32] - v[1163] * v[33];
	v[1270] = -(v[1012] * v[30]) - v[1064] * v[31] - v[1116] * v[32] - v[1168] * v[33];
	v[1271] = -(v[30] * v[582]) - v[31] * v[583] - v[32] * v[584] - v[33] * v[585];
	v[1272] = -(v[30] * v[565]) - v[31] * v[566] - v[32] * v[567] - v[33] * v[568];
	v[1273] = -(v[30] * v[548]) - v[31] * v[549] - v[32] * v[550] - v[33] * v[551];
	v[1274] = -(v[1015] * v[30]) - v[1067] * v[31] - v[1119] * v[32] - v[1171] * v[33];
	v[1275] = -(v[1020] * v[30]) - v[1072] * v[31] - v[1124] * v[32] - v[1176] * v[33];
	v[1276] = -(v[1025] * v[30]) - v[1077] * v[31] - v[1129] * v[32] - v[1181] * v[33];
	v[1277] = -(v[30] * v[586]) - v[31] * v[587] - v[32] * v[589] - v[33] * v[590];
	v[1278] = -(v[30] * v[569]) - v[31] * v[570] - v[32] * v[572] - v[33] * v[573];
	v[1279] = -(v[30] * v[552]) - v[31] * v[553] - v[32] * v[555] - v[33] * v[556];
	v[1280] = -(v[1028] * v[30]) - v[1080] * v[31] - v[1132] * v[32] - v[1184] * v[33];
	v[1281] = -(v[1033] * v[30]) - v[1085] * v[31] - v[1137] * v[32] - v[1189] * v[33];
	v[1282] = -(v[1038] * v[30]) - v[1090] * v[31] - v[1142] * v[32] - v[1194] * v[33];
	v[9902] = v[1259];
	v[9903] = v[1260];
	v[9904] = v[1261];
	v[9905] = v[1262];
	v[9906] = v[1263];
	v[9907] = v[1264];
	v[9908] = v[1265];
	v[9909] = v[1266];
	v[9910] = v[1267];
	v[9911] = v[1268];
	v[9912] = v[1269];
	v[9913] = v[1270];
	v[9914] = v[1271];
	v[9915] = v[1272];
	v[9916] = v[1273];
	v[9917] = v[1274];
	v[9918] = v[1275];
	v[9919] = v[1276];
	v[9920] = v[1277];
	v[9921] = v[1278];
	v[9922] = v[1279];
	v[9923] = v[1280];
	v[9924] = v[1281];
	v[9925] = v[1282];
	v[1728] = v[1259] * v[349] + v[1260] * v[350] + v[1261] * v[351] + v[1265] * v[352] + v[1266] * v[353] + v[1267] * v[354]
		+ v[1271] * v[355] + v[1272] * v[356] + v[1273] * v[357] + v[1277] * v[358] + v[1278] * v[359] + v[1279] * v[360]
		+ v[1262] * v[366] + v[1263] * v[376] + v[1264] * v[386] + v[1268] * v[392] + v[1269] * v[402] + v[1270] * v[412]
		+ v[1274] * v[418] + v[1275] * v[428] + v[1276] * v[438] + v[1280] * v[444] + v[1281] * v[454] + v[1282] * v[464];
	v[1283] = -(v[34] * v[591]) - v[35] * v[592] - v[36] * v[593] - v[37] * v[594];
	v[1568] = v[1211] * v[500] + v[1235] * v[518] - v[1259] * v[528] - v[1283] * v[546];
	v[1567] = v[1211] * v[497] + v[1235] * v[517] - v[1259] * v[525] - v[1283] * v[545];
	v[1566] = v[1211] * v[494] + v[1235] * v[516] - v[1259] * v[522] - v[1283] * v[544];
	v[1284] = -(v[34] * v[574]) - v[35] * v[575] - v[36] * v[576] - v[37] * v[577];
	v[1572] = v[1212] * v[500] + v[1236] * v[518] - v[1260] * v[528] - v[1284] * v[546];
	v[1571] = v[1212] * v[497] + v[1236] * v[517] - v[1260] * v[525] - v[1284] * v[545];
	v[1570] = v[1212] * v[494] + v[1236] * v[516] - v[1260] * v[522] - v[1284] * v[544];
	v[1285] = -(v[34] * v[557]) - v[35] * v[558] - v[36] * v[559] - v[37] * v[560];
	v[1576] = v[1213] * v[500] + v[1237] * v[518] - v[1261] * v[528] - v[1285] * v[546];
	v[1575] = v[1213] * v[497] + v[1237] * v[517] - v[1261] * v[525] - v[1285] * v[545];
	v[1574] = v[1213] * v[494] + v[1237] * v[516] - v[1261] * v[522] - v[1285] * v[544];
	v[1286] = -(v[1041] * v[35]) - v[1093] * v[36] - v[1145] * v[37] - v[34] * v[989];
	v[1580] = v[1214] * v[500] + v[1238] * v[518] - v[1262] * v[528] - v[1286] * v[546];
	v[1579] = v[1214] * v[497] + v[1238] * v[517] - v[1262] * v[525] - v[1286] * v[545];
	v[1578] = v[1214] * v[494] + v[1238] * v[516] - v[1262] * v[522] - v[1286] * v[544];
	v[1287] = -(v[1046] * v[35]) - v[1098] * v[36] - v[1150] * v[37] - v[34] * v[994];
	v[1584] = v[1215] * v[500] + v[1239] * v[518] - v[1263] * v[528] - v[1287] * v[546];
	v[1583] = v[1215] * v[497] + v[1239] * v[517] - v[1263] * v[525] - v[1287] * v[545];
	v[1582] = v[1215] * v[494] + v[1239] * v[516] - v[1263] * v[522] - v[1287] * v[544];
	v[1288] = -(v[1051] * v[35]) - v[1103] * v[36] - v[1155] * v[37] - v[34] * v[999];
	v[1588] = v[1216] * v[500] + v[1240] * v[518] - v[1264] * v[528] - v[1288] * v[546];
	v[1587] = v[1216] * v[497] + v[1240] * v[517] - v[1264] * v[525] - v[1288] * v[545];
	v[1586] = v[1216] * v[494] + v[1240] * v[516] - v[1264] * v[522] - v[1288] * v[544];
	v[1289] = -(v[34] * v[595]) - v[35] * v[596] - v[36] * v[597] - v[37] * v[598];
	v[1592] = v[1217] * v[500] + v[1241] * v[518] - v[1265] * v[528] - v[1289] * v[546];
	v[1591] = v[1217] * v[497] + v[1241] * v[517] - v[1265] * v[525] - v[1289] * v[545];
	v[1590] = v[1217] * v[494] + v[1241] * v[516] - v[1265] * v[522] - v[1289] * v[544];
	v[1290] = -(v[34] * v[578]) - v[35] * v[579] - v[36] * v[580] - v[37] * v[581];
	v[1596] = v[1218] * v[500] + v[1242] * v[518] - v[1266] * v[528] - v[1290] * v[546];
	v[1595] = v[1218] * v[497] + v[1242] * v[517] - v[1266] * v[525] - v[1290] * v[545];
	v[1594] = v[1218] * v[494] + v[1242] * v[516] - v[1266] * v[522] - v[1290] * v[544];
	v[1291] = -(v[34] * v[561]) - v[35] * v[562] - v[36] * v[563] - v[37] * v[564];
	v[1600] = v[1219] * v[500] + v[1243] * v[518] - v[1267] * v[528] - v[1291] * v[546];
	v[1599] = v[1219] * v[497] + v[1243] * v[517] - v[1267] * v[525] - v[1291] * v[545];
	v[1598] = v[1219] * v[494] + v[1243] * v[516] - v[1267] * v[522] - v[1291] * v[544];
	v[1292] = -(v[1002] * v[34]) - v[1054] * v[35] - v[1106] * v[36] - v[1158] * v[37];
	v[1604] = v[1220] * v[500] + v[1244] * v[518] - v[1268] * v[528] - v[1292] * v[546];
	v[1603] = v[1220] * v[497] + v[1244] * v[517] - v[1268] * v[525] - v[1292] * v[545];
	v[1602] = v[1220] * v[494] + v[1244] * v[516] - v[1268] * v[522] - v[1292] * v[544];
	v[1293] = -(v[1007] * v[34]) - v[1059] * v[35] - v[1111] * v[36] - v[1163] * v[37];
	v[1608] = v[1221] * v[500] + v[1245] * v[518] - v[1269] * v[528] - v[1293] * v[546];
	v[1607] = v[1221] * v[497] + v[1245] * v[517] - v[1269] * v[525] - v[1293] * v[545];
	v[1606] = v[1221] * v[494] + v[1245] * v[516] - v[1269] * v[522] - v[1293] * v[544];
	v[1294] = -(v[1012] * v[34]) - v[1064] * v[35] - v[1116] * v[36] - v[1168] * v[37];
	v[1612] = v[1222] * v[500] + v[1246] * v[518] - v[1270] * v[528] - v[1294] * v[546];
	v[1611] = v[1222] * v[497] + v[1246] * v[517] - v[1270] * v[525] - v[1294] * v[545];
	v[1610] = v[1222] * v[494] + v[1246] * v[516] - v[1270] * v[522] - v[1294] * v[544];
	v[1295] = -(v[34] * v[582]) - v[35] * v[583] - v[36] * v[584] - v[37] * v[585];
	v[1616] = v[1223] * v[500] + v[1247] * v[518] - v[1271] * v[528] - v[1295] * v[546];
	v[1615] = v[1223] * v[497] + v[1247] * v[517] - v[1271] * v[525] - v[1295] * v[545];
	v[1614] = v[1223] * v[494] + v[1247] * v[516] - v[1271] * v[522] - v[1295] * v[544];
	v[1296] = -(v[34] * v[565]) - v[35] * v[566] - v[36] * v[567] - v[37] * v[568];
	v[1620] = v[1224] * v[500] + v[1248] * v[518] - v[1272] * v[528] - v[1296] * v[546];
	v[1619] = v[1224] * v[497] + v[1248] * v[517] - v[1272] * v[525] - v[1296] * v[545];
	v[1618] = v[1224] * v[494] + v[1248] * v[516] - v[1272] * v[522] - v[1296] * v[544];
	v[1297] = -(v[34] * v[548]) - v[35] * v[549] - v[36] * v[550] - v[37] * v[551];
	v[1624] = v[1225] * v[500] + v[1249] * v[518] - v[1273] * v[528] - v[1297] * v[546];
	v[1623] = v[1225] * v[497] + v[1249] * v[517] - v[1273] * v[525] - v[1297] * v[545];
	v[1622] = v[1225] * v[494] + v[1249] * v[516] - v[1273] * v[522] - v[1297] * v[544];
	v[1298] = -(v[1015] * v[34]) - v[1067] * v[35] - v[1119] * v[36] - v[1171] * v[37];
	v[1628] = v[1226] * v[500] + v[1250] * v[518] - v[1274] * v[528] - v[1298] * v[546];
	v[1627] = v[1226] * v[497] + v[1250] * v[517] - v[1274] * v[525] - v[1298] * v[545];
	v[1626] = v[1226] * v[494] + v[1250] * v[516] - v[1274] * v[522] - v[1298] * v[544];
	v[1299] = -(v[1020] * v[34]) - v[1072] * v[35] - v[1124] * v[36] - v[1176] * v[37];
	v[1632] = v[1227] * v[500] + v[1251] * v[518] - v[1275] * v[528] - v[1299] * v[546];
	v[1631] = v[1227] * v[497] + v[1251] * v[517] - v[1275] * v[525] - v[1299] * v[545];
	v[1630] = v[1227] * v[494] + v[1251] * v[516] - v[1275] * v[522] - v[1299] * v[544];
	v[1300] = -(v[1025] * v[34]) - v[1077] * v[35] - v[1129] * v[36] - v[1181] * v[37];
	v[1636] = v[1228] * v[500] + v[1252] * v[518] - v[1276] * v[528] - v[1300] * v[546];
	v[1635] = v[1228] * v[497] + v[1252] * v[517] - v[1276] * v[525] - v[1300] * v[545];
	v[1634] = v[1228] * v[494] + v[1252] * v[516] - v[1276] * v[522] - v[1300] * v[544];
	v[1301] = -(v[34] * v[586]) - v[35] * v[587] - v[36] * v[589] - v[37] * v[590];
	v[1640] = v[1229] * v[500] + v[1253] * v[518] - v[1277] * v[528] - v[1301] * v[546];
	v[1639] = v[1229] * v[497] + v[1253] * v[517] - v[1277] * v[525] - v[1301] * v[545];
	v[1638] = v[1229] * v[494] + v[1253] * v[516] - v[1277] * v[522] - v[1301] * v[544];
	v[1302] = -(v[34] * v[569]) - v[35] * v[570] - v[36] * v[572] - v[37] * v[573];
	v[1644] = v[1230] * v[500] + v[1254] * v[518] - v[1278] * v[528] - v[1302] * v[546];
	v[1643] = v[1230] * v[497] + v[1254] * v[517] - v[1278] * v[525] - v[1302] * v[545];
	v[1642] = v[1230] * v[494] + v[1254] * v[516] - v[1278] * v[522] - v[1302] * v[544];
	v[1303] = -(v[34] * v[552]) - v[35] * v[553] - v[36] * v[555] - v[37] * v[556];
	v[1648] = v[1231] * v[500] + v[1255] * v[518] - v[1279] * v[528] - v[1303] * v[546];
	v[1647] = v[1231] * v[497] + v[1255] * v[517] - v[1279] * v[525] - v[1303] * v[545];
	v[1646] = v[1231] * v[494] + v[1255] * v[516] - v[1279] * v[522] - v[1303] * v[544];
	v[1304] = -(v[1028] * v[34]) - v[1080] * v[35] - v[1132] * v[36] - v[1184] * v[37];
	v[1652] = v[1232] * v[500] + v[1256] * v[518] - v[1280] * v[528] - v[1304] * v[546];
	v[1651] = v[1232] * v[497] + v[1256] * v[517] - v[1280] * v[525] - v[1304] * v[545];
	v[1650] = v[1232] * v[494] + v[1256] * v[516] - v[1280] * v[522] - v[1304] * v[544];
	v[1305] = -(v[1033] * v[34]) - v[1085] * v[35] - v[1137] * v[36] - v[1189] * v[37];
	v[1656] = v[1233] * v[500] + v[1257] * v[518] - v[1281] * v[528] - v[1305] * v[546];
	v[1655] = v[1233] * v[497] + v[1257] * v[517] - v[1281] * v[525] - v[1305] * v[545];
	v[1654] = v[1233] * v[494] + v[1257] * v[516] - v[1281] * v[522] - v[1305] * v[544];
	v[1306] = -(v[1038] * v[34]) - v[1090] * v[35] - v[1142] * v[36] - v[1194] * v[37];
	v[9854] = v[1283];
	v[9855] = v[1284];
	v[9856] = v[1285];
	v[9857] = v[1286];
	v[9858] = v[1287];
	v[9859] = v[1288];
	v[9860] = v[1289];
	v[9861] = v[1290];
	v[9862] = v[1291];
	v[9863] = v[1292];
	v[9864] = v[1293];
	v[9865] = v[1294];
	v[9866] = v[1295];
	v[9867] = v[1296];
	v[9868] = v[1297];
	v[9869] = v[1298];
	v[9870] = v[1299];
	v[9871] = v[1300];
	v[9872] = v[1301];
	v[9873] = v[1302];
	v[9874] = v[1303];
	v[9875] = v[1304];
	v[9876] = v[1305];
	v[9877] = v[1306];
	v[1730] = v[1283] * v[349] + v[1284] * v[350] + v[1285] * v[351] + v[1289] * v[352] + v[1290] * v[353] + v[1291] * v[354]
		+ v[1295] * v[355] + v[1296] * v[356] + v[1297] * v[357] + v[1301] * v[358] + v[1302] * v[359] + v[1303] * v[360]
		+ v[1286] * v[366] + v[1287] * v[376] + v[1288] * v[386] + v[1292] * v[392] + v[1293] * v[402] + v[1294] * v[412]
		+ v[1298] * v[418] + v[1299] * v[428] + v[1300] * v[438] + v[1304] * v[444] + v[1305] * v[454] + v[1306] * v[464];
	v[1660] = v[1234] * v[500] + v[1258] * v[518] - v[1282] * v[528] - v[1306] * v[546];
	v[1659] = v[1234] * v[497] + v[1258] * v[517] - v[1282] * v[525] - v[1306] * v[545];
	v[1658] = v[1234] * v[494] + v[1258] * v[516] - v[1282] * v[522] - v[1306] * v[544];
	b1307 = sqrt(Power(v[481] * v[489] - v[480] * v[490], 2) + Power(-(v[482] * v[489]) + v[480] * v[491], 2) + Power
	(v[482] * v[490] - v[481] * v[491], 2)) > 0.1e-7;
	if (b1307) {
		v[1309] = v[482] * v[490] - v[481] * v[491];
		v[1310] = -(v[482] * v[489]) + v[480] * v[491];
		v[1311] = v[481] * v[489] - v[480] * v[490];
		v[1312] = sqrt((v[1309] * v[1309]) + (v[1310] * v[1310]) + (v[1311] * v[1311]));
		v[2690] = 1e0 / (v[1312] * v[1312]);
		v[1876] = v[1312];
		v[2701] = 1e0 - (v[1876] * v[1876]);
		v[5357] = 1e0 / Power(v[2701], 0.15e1);
		v[2696] = 1e0 / sqrt(v[2701]);
		v[1875] = asin(v[1876]) / 2e0;
		v[2695] = 1e0 / Power(cos(v[1875]), 2);
		v[5007] = v[2695] * v[2696];
		v[1314] = 2e0 * tan(v[1875]);
		if (v[1312] > 0.1e-7) { v07 = 1e0 / v[1312]; v08 = (-(v07 / v[1312])); v09 = (2e0 * v07) / (v[1312] * v[1312]); }
		else {
			v07 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[1312]) * (0.7199999994e10 - 0.7199999982e18 * v[1312] + 0.6e25 * Power
			(v[1312], 3) + 0.23999999819999998e26 * (v[1312] * v[1312])));
			v08 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[1312] + 0.6e25 * Power(v[1312], 3)
				+ 0.1799999982e26 * (v[1312] * v[1312]));
			v09 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[1312] - 0.3e17 * (v[1312] * v[1312]));
		};
		v[1318] = v09;
		v[1319] = v08;
		v[1320] = v07;
		v[5356] = v[1314] * v[1319] + v[1320] * v[5007];
		v[4818] = v[1314] * v[1320];
		v[1321] = v[1309] * v[4818];
		v[5042] = 2e0 * v[1321];
		v[4903] = v[1321] / 2e0;
		v[1332] = (v[1321] * v[1321]);
		v[1322] = v[1310] * v[4818];
		v[4819] = v[1322] / 2e0;
		v[1330] = v[1321] * v[4819];
		v[1325] = (v[1322] * v[1322]);
		v[1830] = -v[1325] - v[1332];
		v[1323] = v[1311] * v[4818];
		v[5039] = 2e0 * v[1323];
		v[1860] = -v[1323] + v[1330];
		v[1851] = v[1323] + v[1330];
		v[1337] = v[1323] * v[4819];
		v[1842] = -v[1321] + v[1337];
		v[1834] = v[1321] + v[1337];
		v[1335] = v[1323] * v[4903];
		v[1855] = v[1322] + v[1335];
		v[1838] = -v[1322] + v[1335];
		v[1326] = (v[1323] * v[1323]);
		v[1870] = 4e0 + v[1325] + v[1326] + v[1332];
		v[5358] = 1e0 / Power(v[1870], 3);
		v[5038] = -4e0 / (v[1870] * v[1870]);
		v[1865] = -v[1325] - v[1326];
		v[1847] = -v[1326] - v[1332];
		v[1324] = 4e0 / v[1870];
		v[4820] = v[1324] / 2e0;
		v[1327] = 1e0 + v[1865] * v[4820];
		v[1328] = v[1324] * v[1860];
		v[1329] = v[1324] * v[1855];
		v[1331] = v[1324] * v[1851];
		v[1333] = 1e0 + v[1847] * v[4820];
		v[1334] = v[1324] * v[1842];
		v[1336] = v[1324] * v[1838];
		v[1338] = v[1324] * v[1834];
		v[1339] = 1e0 + v[1830] * v[4820];
	}
	else {
		v[1327] = 1e0;
		v[1328] = 0e0;
		v[1329] = 0e0;
		v[1331] = 0e0;
		v[1333] = 1e0;
		v[1334] = 0e0;
		v[1336] = 0e0;
		v[1338] = 0e0;
		v[1339] = 1e0;
	};
	if (b39) {
		v[1787] = 1e0 - v[1455];
		v[1785] = 1e0 - v[1413];
		v[1783] = 1e0 - v[1383];
		v[1344] = v[10] * v[1338] + v[11] * v[1339] + v[237] * (v[233] + v[219] * v[244] + v[220] * v[246]) + v[238] * (v[236]
			+ v[228] * v[244] + v[229] * v[246]) - v[324] * (v[320] + v[306] * v[331] + v[307] * v[333]) - v[325] * (v[323]
				+ v[315] * v[331] + v[316] * v[333]) + v[1336] * v[9];
		v[4821] = v[1344] * v[482];
		v[1343] = v[10] * v[1333] + v[11] * v[1334] + v[237] * (v[232] + v[216] * v[244] + v[217] * v[246]) + v[238] * (v[235]
			+ v[225] * v[244] + v[226] * v[246]) - v[324] * (v[319] + v[303] * v[331] + v[304] * v[333]) - v[325] * (v[322]
				+ v[312] * v[331] + v[313] * v[333]) + v[1331] * v[9];
		v[4823] = v[1343] * v[481];
		v[4902] = v[4821] + v[4823];
		v[1342] = v[10] * v[1328] + v[11] * v[1329] + v[237] * (v[231] + v[213] * v[244] + v[214] * v[246]) + v[238] * (v[234]
			+ v[222] * v[244] + v[223] * v[246]) - v[324] * (v[318] + v[300] * v[331] + v[301] * v[333]) - v[325] * (v[321]
				+ v[309] * v[331] + v[310] * v[333]) + v[1327] * v[9];
		v[4822] = -(v[1342] * v[480]);
		v[4901] = -v[4821] + v[4822];
		v[4900] = v[4822] - v[4823];
		v[1341] = v[1342] * v[1783] - v[480] * v[4902];
		v[1345] = v[1343] * v[1785] + v[481] * v[4901];
		v[1346] = v[1344] * v[1787] + v[482] * v[4900];
	}
	else {
		v[1341] = 0e0;
		v[1345] = 0e0;
		v[1346] = 0e0;
	};
	v[1384] = -0.5e0 * v[1361] + v[4825];
	v[1421] = v[1359] + v[1360] * v[964] + v[1384] * v[967];
	v[3154] = v[1421] * v[480];
	v[1418] = -v[1360] + v[1384] * v[963] + v[1359] * v[964];
	v[3153] = v[1418] * v[480];
	v[1385] = -0.5e0 * v[1360] + v[4826];
	v[1425] = -v[1359] + v[1361] * v[966] + v[1385] * v[968];
	v[3256] = v[1425] * v[481];
	v[1419] = v[1361] + v[1385] * v[963] + v[1359] * v[966];
	v[3252] = v[1419] * v[481];
	v[4849] = v[3153] + v[3252];
	v[1386] = v[1360] * v[966];
	v[1387] = v[1361] * v[964];
	v[1420] = v[1386] + v[1387] + v[4827] * v[963];
	v[4824] = v[1420] * v[482];
	v[4880] = v[3153] + v[4824];
	v[4868] = v[3252] + v[4824];
	v[3277] = v[1383] * v[1418] + v[480] * v[4868];
	v[3167] = v[1413] * v[1419] + v[481] * v[4880];
	v[3065] = v[1420] * v[1455] + v[482] * v[4849];
	v[1388] = v[1359] * v[965];
	v[1424] = v[1386] + v[1388] + v[4825] * v[968];
	v[3155] = v[1424] * v[480];
	v[4851] = v[3155] + v[3256];
	v[1422] = v[1387] + v[1388] + v[4826] * v[967];
	v[3254] = v[1422] * v[481];
	v[4850] = v[3154] + v[3254];
	v[1389] = -0.5e0 * v[1359] + v[4827];
	v[1426] = v[1360] + v[1361] * v[965] + v[1389] * v[968];
	v[4828] = v[1426] * v[482];
	v[4882] = v[3155] + v[4828];
	v[4870] = v[3256] + v[4828];
	v[3281] = v[1383] * v[1424] + v[480] * v[4870];
	v[3171] = v[1413] * v[1425] + v[481] * v[4882];
	v[3069] = v[1426] * v[1455] + v[482] * v[4851];
	v[1423] = -v[1361] + v[1360] * v[965] + v[1389] * v[967];
	v[4829] = v[1423] * v[482];
	v[4881] = v[3154] + v[4829];
	v[4869] = v[3254] + v[4829];
	v[3279] = v[1383] * v[1421] + v[480] * v[4869];
	v[3169] = v[1413] * v[1422] + v[481] * v[4881];
	v[3067] = v[1423] * v[1455] + v[482] * v[4850];
	v[1390] = -0.5e0 * v[1367] + v[4831];
	v[1430] = v[1365] + v[1366] * v[970] + v[1390] * v[973];
	v[3157] = v[1430] * v[480];
	v[1427] = -v[1366] + v[1390] * v[969] + v[1365] * v[970];
	v[3156] = v[1427] * v[480];
	v[1391] = -0.5e0 * v[1366] + v[4832];
	v[1434] = -v[1365] + v[1367] * v[972] + v[1391] * v[974];
	v[3262] = v[1434] * v[481];
	v[1428] = v[1367] + v[1391] * v[969] + v[1365] * v[972];
	v[3258] = v[1428] * v[481];
	v[4852] = v[3156] + v[3258];
	v[1392] = v[1366] * v[972];
	v[1393] = v[1367] * v[970];
	v[1429] = v[1392] + v[1393] + v[4833] * v[969];
	v[4830] = v[1429] * v[482];
	v[4883] = v[3156] + v[4830];
	v[4871] = v[3258] + v[4830];
	v[3283] = v[1383] * v[1427] + v[480] * v[4871];
	v[3173] = v[1413] * v[1428] + v[481] * v[4883];
	v[3071] = v[1429] * v[1455] + v[482] * v[4852];
	v[1394] = v[1365] * v[971];
	v[1433] = v[1392] + v[1394] + v[4831] * v[974];
	v[3158] = v[1433] * v[480];
	v[4854] = v[3158] + v[3262];
	v[1431] = v[1393] + v[1394] + v[4832] * v[973];
	v[3260] = v[1431] * v[481];
	v[4853] = v[3157] + v[3260];
	v[1395] = -0.5e0 * v[1365] + v[4833];
	v[1435] = v[1366] + v[1367] * v[971] + v[1395] * v[974];
	v[4834] = v[1435] * v[482];
	v[4885] = v[3158] + v[4834];
	v[4873] = v[3262] + v[4834];
	v[3287] = v[1383] * v[1433] + v[480] * v[4873];
	v[3177] = v[1413] * v[1434] + v[481] * v[4885];
	v[3075] = v[1435] * v[1455] + v[482] * v[4854];
	v[1432] = -v[1367] + v[1366] * v[971] + v[1395] * v[973];
	v[4835] = v[1432] * v[482];
	v[4884] = v[3157] + v[4835];
	v[4872] = v[3260] + v[4835];
	v[3285] = v[1383] * v[1430] + v[480] * v[4872];
	v[3175] = v[1413] * v[1431] + v[481] * v[4884];
	v[3073] = v[1432] * v[1455] + v[482] * v[4853];
	v[1396] = -0.5e0 * v[1373] + v[4837];
	v[1439] = v[1371] + v[1372] * v[976] + v[1396] * v[979];
	v[3160] = v[1439] * v[480];
	v[1436] = -v[1372] + v[1396] * v[975] + v[1371] * v[976];
	v[3159] = v[1436] * v[480];
	v[1397] = -0.5e0 * v[1372] + v[4838];
	v[1443] = -v[1371] + v[1373] * v[978] + v[1397] * v[980];
	v[3268] = v[1443] * v[481];
	v[1437] = v[1373] + v[1397] * v[975] + v[1371] * v[978];
	v[3264] = v[1437] * v[481];
	v[4855] = v[3159] + v[3264];
	v[1398] = v[1372] * v[978];
	v[1399] = v[1373] * v[976];
	v[1438] = v[1398] + v[1399] + v[4839] * v[975];
	v[4836] = v[1438] * v[482];
	v[4886] = v[3159] + v[4836];
	v[4874] = v[3264] + v[4836];
	v[3289] = v[1383] * v[1436] + v[480] * v[4874];
	v[3179] = v[1413] * v[1437] + v[481] * v[4886];
	v[3077] = v[1438] * v[1455] + v[482] * v[4855];
	v[1400] = v[1371] * v[977];
	v[1442] = v[1398] + v[1400] + v[4837] * v[980];
	v[3161] = v[1442] * v[480];
	v[4857] = v[3161] + v[3268];
	v[1440] = v[1399] + v[1400] + v[4838] * v[979];
	v[3266] = v[1440] * v[481];
	v[4856] = v[3160] + v[3266];
	v[1401] = -0.5e0 * v[1371] + v[4839];
	v[1444] = v[1372] + v[1373] * v[977] + v[1401] * v[980];
	v[4840] = v[1444] * v[482];
	v[4888] = v[3161] + v[4840];
	v[4876] = v[3268] + v[4840];
	v[3293] = v[1383] * v[1442] + v[480] * v[4876];
	v[3183] = v[1413] * v[1443] + v[481] * v[4888];
	v[3081] = v[1444] * v[1455] + v[482] * v[4857];
	v[1441] = -v[1373] + v[1372] * v[977] + v[1401] * v[979];
	v[4841] = v[1441] * v[482];
	v[4887] = v[3160] + v[4841];
	v[4875] = v[3266] + v[4841];
	v[3291] = v[1383] * v[1439] + v[480] * v[4875];
	v[3181] = v[1413] * v[1440] + v[481] * v[4887];
	v[3079] = v[1441] * v[1455] + v[482] * v[4856];
	v[1402] = -0.5e0 * v[1379] + v[4843];
	v[1448] = v[1377] + v[1378] * v[982] + v[1402] * v[985];
	v[3163] = v[1448] * v[480];
	v[1445] = -v[1378] + v[1402] * v[981] + v[1377] * v[982];
	v[3162] = v[1445] * v[480];
	v[1403] = -0.5e0 * v[1378] + v[4846];
	v[1452] = -v[1377] + v[1379] * v[984] + v[1403] * v[986];
	v[3274] = v[1452] * v[481];
	v[1446] = v[1379] + v[1403] * v[981] + v[1377] * v[984];
	v[3270] = v[1446] * v[481];
	v[4858] = v[3162] + v[3270];
	v[1404] = v[1378] * v[984];
	v[1405] = v[1379] * v[982];
	v[1447] = v[1404] + v[1405] + v[4859] * v[981];
	v[4842] = v[1447] * v[482];
	v[4889] = v[3162] + v[4842];
	v[4877] = v[3270] + v[4842];
	v[3295] = v[1383] * v[1445] + v[480] * v[4877];
	v[3185] = v[1413] * v[1446] + v[481] * v[4889];
	v[3083] = v[1447] * v[1455] + v[482] * v[4858];
	v[1406] = v[1377] * v[983];
	v[1451] = v[1404] + v[1406] + v[4843] * v[986];
	v[4844] = v[1451] * v[366] + v[1448] * v[376] + v[1445] * v[386] + v[1442] * v[392] + v[1439] * v[402] + v[1436] * v[412]
		+ v[1433] * v[418] + v[1430] * v[428] + v[1427] * v[438] + v[1424] * v[444] + v[1421] * v[454] + v[1418] * v[464];
	v[4845] = v[239] * v[349] + v[240] * v[352] - v[326] * v[355] - v[327] * v[358] + v[4844];
	v[4401] = v[482] * v[4845];
	v[4400] = v[481] * v[4844];
	v[3164] = v[1451] * v[480];
	v[4861] = v[3164] + v[3274];
	v[1449] = v[1405] + v[1406] + v[4846] * v[985];
	v[4848] = v[1452] * v[366] + v[1449] * v[376] + v[1446] * v[386] + v[1443] * v[392] + v[1440] * v[402] + v[1437] * v[412]
		+ v[1434] * v[418] + v[1431] * v[428] + v[1428] * v[438] + v[1425] * v[444] + v[1422] * v[454] + v[1419] * v[464];
	v[4847] = v[239] * v[350] + v[240] * v[353] - v[326] * v[356] - v[327] * v[359] + v[4848];
	v[4395] = v[482] * v[4847];
	v[4391] = v[480] * v[4848];
	v[3272] = v[1449] * v[481];
	v[4892] = v[3163] + v[3272];
	v[4389] = v[3165] + v[17] * (v[464] * v[4849] + v[454] * v[4850] + v[444] * v[4851] + v[438] * v[4852] + v[428] * v[4853]
		+ v[418] * v[4854] + v[412] * v[4855] + v[402] * v[4856] + v[392] * v[4857] + v[386] * v[4858] + v[366] * v[4861]
		+ v[376] * v[4892] + v[4893]);
	v[1407] = -0.5e0 * v[1377] + v[4859];
	v[1453] = v[1378] + v[1379] * v[983] + v[1407] * v[986];
	v[4860] = v[1453] * v[482];
	v[4890] = v[3164] + v[4860];
	v[4879] = v[3274] + v[4860];
	v[3299] = v[1383] * v[1451] + v[480] * v[4879];
	v[3189] = v[1413] * v[1452] + v[481] * v[4890];
	v[3087] = v[1453] * v[1455] + v[482] * v[4861];
	v[1450] = -v[1379] + v[1378] * v[983] + v[1407] * v[985];
	v[4867] = v[1450] * v[482];
	v[4878] = v[3272] + v[4867];
	v[4866] = v[1453] * v[366] + v[1450] * v[376] + v[1447] * v[386] + v[1444] * v[392] + v[1441] * v[402] + v[1438] * v[412]
		+ v[1435] * v[418] + v[1432] * v[428] + v[1429] * v[438] + v[1426] * v[444] + v[1423] * v[454] + v[1420] * v[464];
	v[4388] = v[4862] + v[4863] + v[4864] + v[4865] + v[4866];
	v[5046] = v[4388] - v[4866];
	v[4387] = v[481] * v[4866];
	v[4386] = v[480] * v[4866];
	v[3297] = v[1383] * v[1448] + v[480] * v[4878];
	v[4891] = v[3163] + v[4867];
	v[4399] = v[3165] + v[17] * (v[3251] + v[464] * v[4868] + v[454] * v[4869] + v[444] * v[4870] + v[438] * v[4871]
		+ v[428] * v[4872] + v[418] * v[4873] + v[412] * v[4874] + v[402] * v[4875] + v[392] * v[4876] + v[386] * v[4877]
		+ v[376] * v[4878] + v[366] * v[4879]);
	v[4394] = v[3165] + v[17] * (v[3251] + v[464] * v[4880] + v[454] * v[4881] + v[444] * v[4882] + v[438] * v[4883]
		+ v[428] * v[4884] + v[418] * v[4885] + v[412] * v[4886] + v[402] * v[4887] + v[392] * v[4888] + v[386] * v[4889]
		+ v[366] * v[4890] + v[376] * v[4891]);
	v[3187] = v[1413] * v[1449] + v[481] * v[4891];
	v[3085] = v[1450] * v[1455] + v[482] * v[4892];
	v[1483] = v[3165] * v[480] + v[17] * (v[1409] * v[350] + v[1410] * v[353] + v[1411] * v[356] + v[1412] * v[359]
		+ v[3299] * v[366] + v[3297] * v[376] + v[3295] * v[386] + v[3293] * v[392] + v[3291] * v[402] + v[3289] * v[412]
		+ v[3287] * v[418] + v[3285] * v[428] + v[3283] * v[438] + v[3281] * v[444] + v[3279] * v[454] + v[3277] * v[464]
		+ v[3251] * v[480] + v[1383] * (-v[4844] + v[4845]));
	v[5281] = v[1483] * v[240];
	v[5280] = v[1483] * v[239];
	v[5260] = -(v[1483] * v[327]);
	v[5259] = -(v[1483] * v[326]);
	v[1484] = v[3165] * v[481] + v[17] * (v[1409] * v[349] + v[1410] * v[352] + v[1411] * v[355] + v[1412] * v[358]
		+ v[3189] * v[366] + v[3187] * v[376] + v[3185] * v[386] + v[3183] * v[392] + v[3181] * v[402] + v[3179] * v[412]
		+ v[3177] * v[418] + v[3175] * v[428] + v[3173] * v[438] + v[3171] * v[444] + v[3169] * v[454] + v[3167] * v[464]
		+ v[3251] * v[481] + v[1413] * (v[4847] - v[4848]));
	v[5283] = v[1484] * v[240];
	v[5282] = v[1484] * v[239];
	v[5262] = -(v[1484] * v[327]);
	v[5261] = -(v[1484] * v[326]);
	v[1485] = v[3165] * v[482] + v[17] * (v[3087] * v[366] + v[3085] * v[376] + v[3083] * v[386] + v[3081] * v[392]
		+ v[3079] * v[402] + v[3077] * v[412] + v[3075] * v[418] + v[3073] * v[428] + v[3071] * v[438] + v[3069] * v[444]
		+ v[3067] * v[454] + v[3065] * v[464] + v[482] * v[4893] + v[1455] * v[5046]);
	v[5476] = v[1483] * v[223] + v[1484] * v[226] + v[1485] * v[229];
	v[5475] = v[1483] * v[222] + v[1484] * v[225] + v[1485] * v[228];
	v[10108] = 0e0;
	v[10109] = 0e0;
	v[10110] = 0e0;
	v[10111] = 0e0;
	v[10112] = 0e0;
	v[10113] = 0e0;
	v[10114] = v[1483];
	v[10115] = v[1484];
	v[10116] = v[1485];
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
	v[10131] = 0e0;
	v[5473] = v[1483] * v[214] + v[1484] * v[217] + v[1485] * v[220];
	v[5472] = v[1483] * v[213] + v[1484] * v[216] + v[1485] * v[219];
	v[10132] = v[1483];
	v[10133] = v[1484];
	v[10134] = v[1485];
	v[10135] = 0e0;
	v[10136] = 0e0;
	v[10137] = 0e0;
	v[10138] = 0e0;
	v[10139] = 0e0;
	v[10140] = 0e0;
	v[10141] = 0e0;
	v[10142] = 0e0;
	v[10143] = 0e0;
	v[10144] = 0e0;
	v[10145] = 0e0;
	v[10146] = 0e0;
	v[10147] = 0e0;
	v[10148] = 0e0;
	v[10149] = 0e0;
	v[10150] = 0e0;
	v[10151] = 0e0;
	v[10152] = 0e0;
	v[10153] = 0e0;
	v[10154] = 0e0;
	v[10155] = 0e0;
	v[5470] = -(v[1483] * v[310]) - v[1484] * v[313] - v[1485] * v[316];
	v[5469] = -(v[1483] * v[309]) - v[1484] * v[312] - v[1485] * v[315];
	v[10060] = 0e0;
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
	v[10076] = 0e0;
	v[10077] = 0e0;
	v[10078] = -v[1483];
	v[10079] = -v[1484];
	v[10080] = -v[1485];
	v[10081] = 0e0;
	v[10082] = 0e0;
	v[10083] = 0e0;
	v[5467] = -(v[1483] * v[301]) - v[1484] * v[304] - v[1485] * v[307];
	v[5466] = -(v[1483] * v[300]) - v[1484] * v[303] - v[1485] * v[306];
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
	v[10094] = 0e0;
	v[10095] = 0e0;
	v[10096] = -v[1483];
	v[10097] = -v[1484];
	v[10098] = -v[1485];
	v[10099] = 0e0;
	v[10100] = 0e0;
	v[10101] = 0e0;
	v[10102] = 0e0;
	v[10103] = 0e0;
	v[10104] = 0e0;
	v[10105] = 0e0;
	v[10106] = 0e0;
	v[10107] = 0e0;
	v[5285] = v[1485] * v[240];
	v[5284] = v[1485] * v[239];
	v[5264] = -(v[1485] * v[327]);
	v[5263] = -(v[1485] * v[326]);
	v[2854] = (v[1483] * v[1483]) + (v[1484] * v[1484]) + (v[1485] * v[1485]);
	v[1486] = v[1341] * v[14];
	v[1487] = v[1345] * v[14];
	v[1488] = v[1346] * v[14];
	v[1492] = v[1486] - v[18] * (v[1566] * v[349] + v[1570] * v[350] + v[1574] * v[351] + v[1590] * v[352] + v[1594] * v[353]
		+ v[1598] * v[354] + v[1614] * v[355] + v[1618] * v[356] + v[1622] * v[357] + v[1638] * v[358] + v[1642] * v[359]
		+ v[1646] * v[360] + v[1578] * v[366] + v[1582] * v[376] + v[1586] * v[386] + v[1602] * v[392] + v[1606] * v[402]
		+ v[1610] * v[412] + v[1626] * v[418] + v[1630] * v[428] + v[1634] * v[438] + v[1650] * v[444] + v[1654] * v[454]
		+ v[1658] * v[464]);
	v[1493] = v[1487] - v[18] * (v[1567] * v[349] + v[1571] * v[350] + v[1575] * v[351] + v[1591] * v[352] + v[1595] * v[353]
		+ v[1599] * v[354] + v[1615] * v[355] + v[1619] * v[356] + v[1623] * v[357] + v[1639] * v[358] + v[1643] * v[359]
		+ v[1647] * v[360] + v[1579] * v[366] + v[1583] * v[376] + v[1587] * v[386] + v[1603] * v[392] + v[1607] * v[402]
		+ v[1611] * v[412] + v[1627] * v[418] + v[1631] * v[428] + v[1635] * v[438] + v[1651] * v[444] + v[1655] * v[454]
		+ v[1659] * v[464]);
	v[1494] = v[1488] - v[18] * (v[1568] * v[349] + v[1572] * v[350] + v[1576] * v[351] + v[1592] * v[352] + v[1596] * v[353]
		+ v[1600] * v[354] + v[1616] * v[355] + v[1620] * v[356] + v[1624] * v[357] + v[1640] * v[358] + v[1644] * v[359]
		+ v[1648] * v[360] + v[1580] * v[366] + v[1584] * v[376] + v[1588] * v[386] + v[1604] * v[392] + v[1608] * v[402]
		+ v[1612] * v[412] + v[1628] * v[418] + v[1632] * v[428] + v[1636] * v[438] + v[1652] * v[444] + v[1656] * v[454]
		+ v[1660] * v[464]);
	v[2850] = (v[1492] * v[1492]) + (v[1493] * v[1493]) + (v[1494] * v[1494]);
	if (b38) {
		b1496 = sqrt((v[1492] * v[1492]) + (v[1493] * v[1493]) + (v[1494] * v[1494])) <= v[15] * sqrt((v[1483] * v[1483]) +
			(v[1484] * v[1484]) + (v[1485] * v[1485]));
		if (b1496) {
			v[1498] = v[1492];
			v[1499] = v[1493];
			v[1500] = v[1494];
			v[1501] = 1e0;
		}
		else {
			v[4894] = v[16] * sqrt(v[2854]);
			v[1502] = sqrt(v[2850]);
			if (v[1502] > 0.1e-5) {
				v010 = 1e0 / v[1502]; v011 = (-(v010 / v[1502])); v012 = (2e0 * v010) / (v[1502] * v[1502]
					);
			}
			else {
				v010 = (24000000e0 - (-1e0 + 1000000e0 * v[1502]) * (71999994e0 - 0.71999982e14 * v[1502] + 0.6e19 * Power(v[1502]
					, 3) + 0.23999982e20 * (v[1502] * v[1502]))) / 24e0;
				v011 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[1502] + 0.6e19 * Power(v[1502], 3) + 0.17999982e20 *
					(v[1502] * v[1502]));
				v012 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[1502] - 0.3e13 * (v[1502] * v[1502]));
			};
			v[1506] = v011;
			v[1507] = v010;
			v[1508] = v[1492] * v[1507];
			v[1509] = v[1493] * v[1507];
			v[1510] = v[1494] * v[1507];
			v[1498] = v[1508] * v[4894];
			v[1499] = v[1509] * v[4894];
			v[1500] = v[1510] * v[4894];
			v[1501] = 0e0;
		};
		if (sqrt((v[1486] * v[1486]) + (v[1487] * v[1487]) + (v[1488] * v[1488])) > v[15] * sqrt((v[1483] * v[1483]) +
			(v[1484] * v[1484]) + (v[1485] * v[1485]))) {
			if (v[14] > 0.1e-5) { v013 = 1e0 / v[14]; v014 = (-(v013 / v[14])); v015 = (2e0 * v013) / (v[14] * v[14]); }
			else {
				v013 = (24000000e0 - (-1e0 + 1000000e0 * v[14]) * (71999994e0 - 0.71999982e14 * v[14] + 0.6e19 * Power(v[14], 3)
					+ 0.23999982e20 * (v[14] * v[14]))) / 24e0;
				v014 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[14] + 0.6e19 * Power(v[14], 3) + 0.17999982e20 *
					(v[14] * v[14]));
				v015 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[14] - 0.3e13 * (v[14] * v[14]));
			};
			v[1520] = sqrt((v[1486] * v[1486]) + (v[1487] * v[1487]) + (v[1488] * v[1488]));
			if (v[1520] > 0.1e-5) {
				v016 = 1e0 / v[1520]; v017 = (-(v016 / v[1520])); v018 = (2e0 * v016) / (v[1520] * v[1520]
					);
			}
			else {
				v016 = (24000000e0 - (-1e0 + 1000000e0 * v[1520]) * (71999994e0 - 0.71999982e14 * v[1520] + 0.6e19 * Power(v[1520]
					, 3) + 0.23999982e20 * (v[1520] * v[1520]))) / 24e0;
				v017 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[1520] + 0.6e19 * Power(v[1520], 3) + 0.17999982e20 *
					(v[1520] * v[1520]));
				v018 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[1520] - 0.3e13 * (v[1520] * v[1520]));
			};
			v[1527] = -(v013 * v016 * v[16] * sqrt(v[2854]));
			v[1526] = v[1341] + v[1486] * v[1527];
			v[1528] = v[1345] + v[1487] * v[1527];
			v[1529] = v[1346] + v[1488] * v[1527];
		}
		else {
			v[1526] = 0e0;
			v[1528] = 0e0;
			v[1529] = 0e0;
		};
	}
	else {
		b1530 = sqrt((v[1492] * v[1492]) + (v[1493] * v[1493]) + (v[1494] * v[1494])) <= v[16] * sqrt((v[1483] * v[1483]) +
			(v[1484] * v[1484]) + (v[1485] * v[1485]));
		if (b1530) {
			v[1498] = v[1492];
			v[1499] = v[1493];
			v[1500] = v[1494];
			v[1501] = 1e0;
		}
		else {
			v[1541] = sqrt(v[2854]);
			v[4895] = v[1541] * v[16];
			v[1532] = sqrt(v[2850]);
			if (v[1532] > 0.1e-5) {
				v019 = 1e0 / v[1532]; v020 = (-(v019 / v[1532])); v021 = (2e0 * v019) / (v[1532] * v[1532]
					);
			}
			else {
				v019 = (24000000e0 - (-1e0 + 1000000e0 * v[1532]) * (71999994e0 - 0.71999982e14 * v[1532] + 0.6e19 * Power(v[1532]
					, 3) + 0.23999982e20 * (v[1532] * v[1532]))) / 24e0;
				v020 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[1532] + 0.6e19 * Power(v[1532], 3) + 0.17999982e20 *
					(v[1532] * v[1532]));
				v021 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[1532] - 0.3e13 * (v[1532] * v[1532]));
			};
			v[1536] = v020;
			v[1537] = v019;
			v[1538] = v[1492] * v[1537];
			v[1539] = v[1493] * v[1537];
			v[1540] = v[1494] * v[1537];
			v[1498] = v[1538] * v[4895];
			v[1499] = v[1539] * v[4895];
			v[1500] = v[1540] * v[4895];
			v[1501] = 0e0;
		};
		if (sqrt((v[1486] * v[1486]) + (v[1487] * v[1487]) + (v[1488] * v[1488])) > v[16] * sqrt((v[1483] * v[1483]) +
			(v[1484] * v[1484]) + (v[1485] * v[1485]))) {
			if (v[14] > 0.1e-5) { v022 = 1e0 / v[14]; v023 = (-(v022 / v[14])); v024 = (2e0 * v022) / (v[14] * v[14]); }
			else {
				v022 = (24000000e0 - (-1e0 + 1000000e0 * v[14]) * (71999994e0 - 0.71999982e14 * v[14] + 0.6e19 * Power(v[14], 3)
					+ 0.23999982e20 * (v[14] * v[14]))) / 24e0;
				v023 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[14] + 0.6e19 * Power(v[14], 3) + 0.17999982e20 *
					(v[14] * v[14]));
				v024 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[14] - 0.3e13 * (v[14] * v[14]));
			};
			v[1550] = sqrt((v[1486] * v[1486]) + (v[1487] * v[1487]) + (v[1488] * v[1488]));
			if (v[1550] > 0.1e-5) {
				v025 = 1e0 / v[1550]; v026 = (-(v025 / v[1550])); v027 = (2e0 * v025) / (v[1550] * v[1550]
					);
			}
			else {
				v025 = (24000000e0 - (-1e0 + 1000000e0 * v[1550]) * (71999994e0 - 0.71999982e14 * v[1550] + 0.6e19 * Power(v[1550]
					, 3) + 0.23999982e20 * (v[1550] * v[1550]))) / 24e0;
				v026 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[1550] + 0.6e19 * Power(v[1550], 3) + 0.17999982e20 *
					(v[1550] * v[1550]));
				v027 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[1550] - 0.3e13 * (v[1550] * v[1550]));
			};
			v[1556] = -(v022 * v025 * v[16] * sqrt(v[2854]));
			v[1526] = v[1341] + v[1486] * v[1556];
			v[1528] = v[1345] + v[1487] * v[1556];
			v[1529] = v[1346] + v[1488] * v[1556];
		}
		else {
			v[1526] = 0e0;
			v[1528] = 0e0;
			v[1529] = 0e0;
		};
	};
	fn[0] = v[1483];
	fn[1] = v[1484];
	fn[2] = v[1485];
	ft[0] = v[1498];
	ft[1] = v[1499];
	ft[2] = v[1500];
	(*stickupdated) = v[1501];
	gtpupdated[0] = v[1341] - v[1526];
	gtpupdated[1] = v[1345] - v[1528];
	gtpupdated[2] = v[1346] - v[1529];
	v[1665] = 0e0;
	v[1666] = 0e0;
	v[1667] = 0e0;
	b1668 = b38;
	if (b1668) {
		b1669 = b1496;
		if (b1669) {
			v[1667] = 0e0;
			v[1666] = 0e0;
			v[1665] = 0e0;
		}
		else {
		};
	}
	else {
		b1670 = b1530;
		if (b1670) {
			v[1667] = 0e0;
			v[1666] = 0e0;
			v[1665] = 0e0;
		}
		else {
		};
	};
	v[4899] = v[1665] * v[18];
	v[2828] = v[1665] * v[4896];
	v[4898] = v[1666] * v[18];
	v[2827] = v[1666] * v[4896];
	v[4897] = v[1667] * v[18];
	v[2826] = v[1667] * v[4896];
	v[2008] = -(v[1728] * v[2826]);
	v[1700] = v[1730] * v[4897];
	v[4908] = v[1700] * v[326];
	v[4907] = v[1700] * v[327];
	v[2105] = -(v[1732] * v[2826]);
	v[1702] = v[1734] * v[4897];
	v[4940] = v[1702] * v[239];
	v[4939] = v[1702] * v[240];
	v[2005] = -(v[1728] * v[2827]);
	v[1731] = v[1730] * v[4898];
	v[4910] = v[1731] * v[326];
	v[4909] = v[1731] * v[327];
	v[2102] = -(v[1732] * v[2827]);
	v[1735] = v[1734] * v[4898];
	v[4942] = v[1735] * v[239];
	v[4941] = v[1735] * v[240];
	v[1737] = -((v[1658] * v[1665] + v[1659] * v[1666] + v[1660] * v[1667]) * v[18]);
	v[2037] = v[1737] * v[2230];
	v[2017] = v[1737] * v[3998];
	v[1738] = -((v[1654] * v[1665] + v[1655] * v[1666] + v[1656] * v[1667]) * v[18]);
	v[2027] = v[1738] * v[3996];
	v[4914] = v[2027] + v[1737] * v[2229];
	v[4913] = v[2017] + v[1738] * v[2227];
	v[1739] = -((v[1650] * v[1665] + v[1651] * v[1666] + v[1652] * v[1667]) * v[18]);
	v[2035] = v[1739] * v[3994];
	v[5220] = v[2035] + v[2037];
	v[4915] = v[2035] + v[1738] * v[2228];
	v[5211] = v[2037] + v[4915];
	v[2029] = v[1739] * v[4765];
	v[5217] = v[2029] + v[4914];
	v[5212] = v[2027] + v[2029];
	v[2020] = v[1739] * v[4767];
	v[5221] = v[2020] + v[4913];
	v[5213] = v[2017] + v[2020];
	v[1740] = -((v[1634] * v[1665] + v[1635] * v[1666] + v[1636] * v[1667]) * v[18]);
	v[2064] = v[1740] * v[2224];
	v[2044] = v[1740] * v[3992];
	v[1741] = -((v[1630] * v[1665] + v[1631] * v[1666] + v[1632] * v[1667]) * v[18]);
	v[2054] = v[1741] * v[3990];
	v[4917] = v[2054] + v[1740] * v[2223];
	v[4916] = v[2044] + v[1741] * v[2221];
	v[1742] = -((v[1626] * v[1665] + v[1627] * v[1666] + v[1628] * v[1667]) * v[18]);
	v[2062] = v[1742] * v[3988];
	v[5207] = v[2062] + v[2064];
	v[4918] = v[2062] + v[1741] * v[2222];
	v[5198] = v[2064] + v[4918];
	v[2056] = v[1742] * v[4760];
	v[5204] = v[2056] + v[4917];
	v[5199] = v[2054] + v[2056];
	v[2047] = v[1742] * v[4762];
	v[5208] = v[2047] + v[4916];
	v[5200] = v[2044] + v[2047];
	v[1743] = -((v[1610] * v[1665] + v[1611] * v[1666] + v[1612] * v[1667]) * v[18]);
	v[2134] = v[1743] * v[2218];
	v[2114] = v[1743] * v[3986];
	v[1744] = -((v[1606] * v[1665] + v[1607] * v[1666] + v[1608] * v[1667]) * v[18]);
	v[2124] = v[1744] * v[3984];
	v[4946] = v[2124] + v[1743] * v[2217];
	v[4945] = v[2114] + v[1744] * v[2215];
	v[1745] = -((v[1602] * v[1665] + v[1603] * v[1666] + v[1604] * v[1667]) * v[18]);
	v[2132] = v[1745] * v[3982];
	v[5194] = v[2132] + v[2134];
	v[4947] = v[2132] + v[1744] * v[2216];
	v[5185] = v[2134] + v[4947];
	v[2126] = v[1745] * v[4755];
	v[5191] = v[2126] + v[4946];
	v[5186] = v[2124] + v[2126];
	v[2117] = v[1745] * v[4757];
	v[5195] = v[2117] + v[4945];
	v[5187] = v[2114] + v[2117];
	v[1746] = -((v[1586] * v[1665] + v[1587] * v[1666] + v[1588] * v[1667]) * v[18]);
	v[2161] = v[1746] * v[2212];
	v[2141] = v[1746] * v[3980];
	v[1747] = -((v[1582] * v[1665] + v[1583] * v[1666] + v[1584] * v[1667]) * v[18]);
	v[2151] = v[1747] * v[3978];
	v[4949] = v[2151] + v[1746] * v[2211];
	v[4948] = v[2141] + v[1747] * v[2209];
	v[1748] = -((v[1578] * v[1665] + v[1579] * v[1666] + v[1580] * v[1667]) * v[18]);
	v[2159] = v[1748] * v[3976];
	v[5181] = v[2159] + v[2161];
	v[4950] = v[2159] + v[1747] * v[2210];
	v[5172] = v[2161] + v[4950];
	v[2153] = v[1748] * v[4750];
	v[5178] = v[2153] + v[4949];
	v[5173] = v[2151] + v[2153];
	v[2144] = v[1748] * v[4752];
	v[5182] = v[2144] + v[4948];
	v[5174] = v[2141] + v[2144];
	v[2002] = -(v[1728] * v[2828]);
	v[1762] = v[1730] * v[4899];
	v[4912] = v[1762] * v[326];
	v[4911] = v[1762] * v[327];
	v[2099] = -(v[1732] * v[2828]);
	v[1764] = v[1734] * v[4899];
	v[4944] = v[1764] * v[239];
	v[4943] = v[1764] * v[240];
	v[1765] = v[14] * v[1667];
	v[2754] = -(v[1765] * v[482]);
	v[1766] = v[14] * v[1666];
	v[2758] = -(v[1766] * v[481]);
	v[2760] = v[2754] + v[2758];
	v[1767] = v[14] * v[1665];
	v[2759] = -(v[1767] * v[480]);
	v[2762] = v[2758] + v[2759];
	v[2761] = v[2754] + v[2759];
	v[1768] = 0e0;
	v[1769] = 0e0;
	v[1770] = 0e0;
	v[1771] = 0e0;
	v[1772] = 0e0;
	v[1773] = 0e0;
	v[1774] = 0e0;
	v[1775] = 0e0;
	v[1776] = 0e0;
	v[1777] = 0e0;
	v[1778] = 0e0;
	v[1779] = 0e0;
	b1780 = b39;
	if (b1780) {
		v[1781] = -(v[1344] * v[1765]);
		v[1782] = -(v[1343] * v[1766]);
		v[1784] = v[1767] * v[1783] + v[2760] * v[480];
		v[1786] = v[1766] * v[1785] + v[2761] * v[481];
		v[1788] = v[1765] * v[1787] + v[2762] * v[482];
		v[1770] = v[1344] * v[2762] + v[1765] * v[4900];
		v[1769] = v[1343] * v[2761] + v[1766] * v[4901];
		v[1790] = -(v[1342] * v[1767]);
		v[1768] = v[1342] * v[2760] - v[1767] * v[4902];
		v[1771] = v[1784] * v[9];
		v[1772] = v[10] * v[1784];
		v[1773] = v[11] * v[1784];
		v[1793] = -(v[1784] * v[325]);
		v[1794] = -(v[1784] * v[324]);
		v[1795] = v[1793] * v[333];
		v[1796] = v[1793] * v[331];
		v[1797] = v[1794] * v[333];
		v[1798] = v[1794] * v[331];
		v[1799] = v[1784] * v[238];
		v[1800] = v[1784] * v[237];
		v[1801] = v[1799] * v[246];
		v[1802] = v[1799] * v[244];
		v[1803] = v[1800] * v[246];
		v[1804] = v[1800] * v[244];
		v[1774] = v[1786] * v[9];
		v[1775] = v[10] * v[1786];
		v[1776] = v[11] * v[1786];
		v[1805] = -(v[1786] * v[325]);
		v[1806] = -(v[1786] * v[324]);
		v[1807] = v[1805] * v[333];
		v[1808] = v[1805] * v[331];
		v[1809] = v[1806] * v[333];
		v[1810] = v[1806] * v[331];
		v[1811] = v[1786] * v[238];
		v[1812] = v[1786] * v[237];
		v[1813] = v[1811] * v[246];
		v[1814] = v[1811] * v[244];
		v[1815] = v[1812] * v[246];
		v[1816] = v[1812] * v[244];
		v[1777] = v[1788] * v[9];
		v[1778] = v[10] * v[1788];
		v[1779] = v[11] * v[1788];
		v[1817] = -(v[1788] * v[325]);
		v[1818] = -(v[1788] * v[324]);
		v[1819] = v[1817] * v[333];
		v[1820] = v[1817] * v[331];
		v[1821] = v[1818] * v[333];
		v[1822] = v[1818] * v[331];
		v[1823] = v[1788] * v[238];
		v[1824] = v[1788] * v[237];
		v[1825] = v[1823] * v[246];
		v[1826] = v[1823] * v[244];
		v[1827] = v[1824] * v[246];
		v[1828] = v[1824] * v[244];
	}
	else {
		v[1804] = 0e0;
		v[1803] = 0e0;
		v[1816] = 0e0;
		v[1815] = 0e0;
		v[1828] = 0e0;
		v[1827] = 0e0;
		v[1802] = 0e0;
		v[1801] = 0e0;
		v[1814] = 0e0;
		v[1813] = 0e0;
		v[1826] = 0e0;
		v[1825] = 0e0;
		v[1800] = 0e0;
		v[1812] = 0e0;
		v[1824] = 0e0;
		v[1799] = 0e0;
		v[1811] = 0e0;
		v[1823] = 0e0;
		v[1798] = 0e0;
		v[1797] = 0e0;
		v[1810] = 0e0;
		v[1809] = 0e0;
		v[1822] = 0e0;
		v[1821] = 0e0;
		v[1796] = 0e0;
		v[1795] = 0e0;
		v[1808] = 0e0;
		v[1807] = 0e0;
		v[1820] = 0e0;
		v[1819] = 0e0;
		v[1794] = 0e0;
		v[1806] = 0e0;
		v[1818] = 0e0;
		v[1793] = 0e0;
		v[1805] = 0e0;
		v[1817] = 0e0;
		v[1790] = 0e0;
		v[1782] = 0e0;
		v[1781] = 0e0;
	};
	v[5004] = 2e0 * v[1790];
	v[5005] = 2e0 * v[1782];
	v[5006] = 2e0 * v[1781];
	v[5013] = v[1771] / 2e0;
	v[5014] = v[1775] / 2e0;
	v[5015] = v[1779] / 2e0;
	b1829 = b1307;
	if (b1829) {
		v[1863] = -(v[1324] * v[1772]);
		v[1858] = v[1324] * v[1773];
		v[1845] = v[1324] * v[1776];
		v[1832] = -(v[1779] * v[4820]);
		v[1836] = v[1324] * v[1778];
		v[1840] = v[1324] * v[1777];
		v[1844] = v[1836] + v[1845];
		v[1849] = -(v[1775] * v[4820]);
		v[1853] = v[1324] * v[1774];
		v[1857] = v[1840] + v[1858];
		v[1864] = v[1853] - v[1863];
		v[1866] = v[1778] * v[1834] + v[1777] * v[1838] + v[1776] * v[1842] + v[1774] * v[1851] + v[1773] * v[1855]
			+ v[1772] * v[1860] + v[1865] * v[5013] + v[1847] * v[5014] + v[1830] * v[5015];
		v[2714] = v[1832] + v[1849] - (4e0 * v[1866]) / (v[1870] * v[1870]);
		v[5012] = 4e0 * v[2714];
		v[2712] = -v[1832] + v[2714] - v[1771] * v[4820];
		v[5011] = 4e0 * (v[1832] - v[1849] + v[2712]);
		v[1871] = v[1853] + v[1863] + 2e0 * v[1323] * v[2712] + v[1844] * v[4819] + v[1857] * v[4903];
		v[1873] = (-2e0 * v[1840] + v[1323] * v[1844] + 2e0 * v[1858] + v[1321] * v[1864] + v[1322] * v[5011]) / 2e0;
		v[1874] = (2e0 * v[1836] - 2e0 * v[1845] + v[1323] * v[1857] + v[1322] * v[1864] + v[1321] * v[5012]) / 2e0;
		v[4904] = v[1311] * v[1871] + v[1310] * v[1873] + v[1309] * v[1874];
		v[2700] = v[1320] * v[4904];
		v[2697] = v[1314] * v[4904];
		v[1877] = v[1319] * v[2697] + v[2700] / (Power(cos(v[1875]), 2) * sqrt(v[2701]));
		v[4905] = v[1877] / v[1312];
		v[1878] = v[1871] * v[4818] + v[1311] * v[4905];
		v[1880] = v[1873] * v[4818] + v[1310] * v[4905];
		v[1881] = v[1874] * v[4818] + v[1309] * v[4905];
		v[1768] = v[1768] - v[1878] * v[490] + v[1880] * v[491];
		v[1770] = v[1770] - v[1880] * v[489] + v[1881] * v[490];
		v[1769] = v[1769] + v[1878] * v[489] - v[1881] * v[491];
	}
	else {
	};
	v[1770] = v[1770] + v[482] * v[5006];
	v[1769] = v[1769] + v[481] * v[5005];
	v[1768] = v[1768] + v[480] * v[5004];
	v[1886] = v[1768] * v[465] + v[1769] * v[466] + v[1770] * v[467];
	v[1888] = v[1886] * v[478];
	v[4906] = v[1888] / v[1475];
	v[1889] = v[1770] * v[479] + v[467] * v[4906];
	v[1891] = v[1769] * v[479] + v[466] * v[4906];
	v[1892] = v[1768] * v[479] + v[465] * v[4906];
	v[1893] = v[1737] * v[2821];
	v[1894] = v[1737] * v[2820];
	v[1895] = v[1737] * v[2819];
	v[1898] = v[1738] * v[2816];
	v[1899] = v[1737] * v[451] + v[1738] * v[453];
	v[4919] = v[1899] * v[286];
	v[1902] = v[1738] * v[2814];
	v[1903] = -v[1893] + v[1898];
	v[1904] = v[1893] + v[1898];
	v[1906] = v[1738] * v[2815] + v[4919] / v[442];
	v[1909] = v[1739] * v[2809];
	v[1910] = v[1737] * v[445] + v[1739] * v[453];
	v[4920] = v[1910] * v[291];
	v[1911] = v[1738] * v[445] + v[1739] * v[451];
	v[4921] = v[1911] * v[295];
	v[5450] = -(v[4764] * v[4919]) - v[441] * v[4920] - v[440] * v[4921] + v[1737] * (-(v[453] * v[455]) - v[445] * v[5114]
		- v[451] * v[5115]) + v[1738] * (-(v[446] * v[451]) - v[445] * v[5116] - v[453] * v[5117]) + v[1739] * (-(v[443] * v[445])
			- v[451] * v[5118] - v[453] * v[5119]);
	v[1913] = v[1739] * v[2811] + v[4920] / v[442];
	v[1914] = v[1906] + v[1913];
	v[1915] = -v[1906] + v[1913];
	v[1917] = v[1739] * v[2810] + v[4921] / v[442];
	v[1918] = v[1894] + v[1917];
	v[1919] = -v[1894] + v[1917];
	v[1920] = v[1740] * v[2806];
	v[1921] = v[1740] * v[2805];
	v[1922] = v[1740] * v[2804];
	v[1925] = v[1741] * v[2801];
	v[1926] = v[1740] * v[425] + v[1741] * v[427];
	v[4929] = v[1926] * v[267];
	v[1929] = v[1741] * v[2799];
	v[1930] = -v[1920] + v[1925];
	v[1931] = v[1920] + v[1925];
	v[1933] = v[1741] * v[2800] + v[4929] / v[416];
	v[1936] = v[1742] * v[2794];
	v[1937] = v[1740] * v[419] + v[1742] * v[427];
	v[4930] = v[1937] * v[272];
	v[1938] = v[1741] * v[419] + v[1742] * v[425];
	v[4931] = v[1938] * v[276];
	v[5453] = -(v[4759] * v[4929]) - v[415] * v[4930] - v[414] * v[4931] + v[1740] * (-(v[427] * v[429]) - v[419] * v[5127]
		- v[425] * v[5128]) + v[1741] * (-(v[420] * v[425]) - v[419] * v[5129] - v[427] * v[5130]) + v[1742] * (-(v[417] * v[419])
			- v[425] * v[5131] - v[427] * v[5132]);
	v[1940] = v[1742] * v[2796] + v[4930] / v[416];
	v[1941] = v[1933] + v[1940];
	v[1942] = -v[1933] + v[1940];
	v[1944] = v[1742] * v[2795] + v[4931] / v[416];
	v[1945] = v[1921] + v[1944];
	v[1946] = -v[1921] + v[1944];
	v[1947] = v[1743] * v[2791];
	v[1948] = v[1743] * v[2790];
	v[1949] = v[1743] * v[2789];
	v[1952] = v[1744] * v[2786];
	v[1953] = v[1743] * v[399] + v[1744] * v[401];
	v[4951] = v[1953] * v[199];
	v[1956] = v[1744] * v[2784];
	v[1957] = -v[1947] + v[1952];
	v[1958] = v[1947] + v[1952];
	v[1960] = v[1744] * v[2785] + v[4951] / v[390];
	v[1963] = v[1745] * v[2779];
	v[1964] = v[1743] * v[393] + v[1745] * v[401];
	v[4952] = v[1964] * v[204];
	v[1965] = v[1744] * v[393] + v[1745] * v[399];
	v[4953] = v[1965] * v[208];
	v[5456] = -(v[4754] * v[4951]) - v[389] * v[4952] - v[388] * v[4953] + v[1743] * (-(v[401] * v[403]) - v[393] * v[5140]
		- v[399] * v[5141]) + v[1744] * (-(v[394] * v[399]) - v[393] * v[5142] - v[401] * v[5143]) + v[1745] * (-(v[391] * v[393])
			- v[399] * v[5144] - v[401] * v[5145]);
	v[1967] = v[1745] * v[2781] + v[4952] / v[390];
	v[1968] = v[1960] + v[1967];
	v[1969] = -v[1960] + v[1967];
	v[1971] = v[1745] * v[2780] + v[4953] / v[390];
	v[1972] = v[1948] + v[1971];
	v[1973] = -v[1948] + v[1971];
	v[1974] = v[1746] * v[2776];
	v[1975] = v[1746] * v[2775];
	v[1976] = v[1746] * v[2774];
	v[1979] = v[1747] * v[2771];
	v[1980] = v[1746] * v[373] + v[1747] * v[375];
	v[4961] = v[180] * v[1980];
	v[1983] = v[1747] * v[2769];
	v[1984] = -v[1974] + v[1979];
	v[1985] = v[1974] + v[1979];
	v[1987] = v[1747] * v[2770] + v[4961] / v[364];
	v[1990] = v[1748] * v[2764];
	v[1991] = v[1746] * v[367] + v[1748] * v[375];
	v[4962] = v[185] * v[1991];
	v[1992] = v[1747] * v[367] + v[1748] * v[373];
	v[4963] = v[189] * v[1992];
	v[5459] = -(v[4749] * v[4961]) - v[363] * v[4962] - v[362] * v[4963] + v[1746] * (-(v[375] * v[377]) - v[367] * v[5153]
		- v[373] * v[5154]) + v[1747] * (-(v[368] * v[373]) - v[367] * v[5155] - v[375] * v[5156]) + v[1748] * (-(v[365] * v[367])
			- v[373] * v[5157] - v[375] * v[5158]);
	v[1994] = v[1748] * v[2766] + v[4962] / v[364];
	v[1995] = v[1987] + v[1994];
	v[1996] = -v[1987] + v[1994];
	v[1998] = v[1748] * v[2765] + v[4963] / v[364];
	v[1999] = v[1975] + v[1998];
	v[2000] = -v[1975] + v[1998];
	v[2001] = v[2002] - v[1892] * v[326];
	v[2003] = -v[2002] - v[1892] * v[327];
	v[1794] = v[1794] + v[2001];
	v[1793] = v[1793] + v[2003];
	v[2004] = v[2005] - v[1891] * v[326];
	v[2006] = -v[2005] - v[1891] * v[327];
	v[1806] = v[1806] + v[2004];
	v[1805] = v[1805] + v[2006];
	v[2007] = v[2008] - v[1889] * v[326];
	v[2009] = -v[2008] - v[1889] * v[327];
	v[1818] = v[1818] + v[2007];
	v[1817] = v[1817] + v[2009];
	v[1819] = v[1819] + v[2009] * v[339] + v[4907] * v[537];
	v[1820] = v[1820] + v[2009] * v[337] + v[4907] * v[536];
	v[1821] = v[1821] + v[2007] * v[339] + v[4908] * v[537];
	v[1822] = v[1822] + v[2007] * v[337] + v[4908] * v[536];
	v[1807] = v[1807] + v[2006] * v[339] + v[4909] * v[537];
	v[1808] = v[1808] + v[2006] * v[337] + v[4909] * v[536];
	v[1809] = v[1809] + v[2004] * v[339] + v[4910] * v[537];
	v[1810] = v[1810] + v[2004] * v[337] + v[4910] * v[536];
	v[1795] = v[1795] + v[2003] * v[339] + v[4911] * v[537];
	v[1796] = v[1796] + v[2003] * v[337] + v[4911] * v[536];
	v[1797] = v[1797] + v[2001] * v[339] + v[4912] * v[537];
	v[1798] = v[1798] + v[2001] * v[337] + v[4912] * v[536];
	v[2022] = v[164] * v[1819] + v[163] * v[1820] + v[453] * v[5221];
	v[2023] = v[161] * v[1819] + v[160] * v[1820] + v[1911] * v[4767] + v[451] * v[4913];
	v[2024] = v[158] * v[1819] + v[157] * v[1820] + v[445] * v[5213];
	v[2031] = v[164] * v[1807] + v[163] * v[1808] + v[1910] * v[4765] + v[453] * v[4914];
	v[2032] = v[161] * v[1807] + v[160] * v[1808] + v[451] * v[5217];
	v[2033] = v[158] * v[1807] + v[157] * v[1808] + v[445] * v[5212];
	v[2040] = v[164] * v[1795] + v[163] * v[1796] + v[1899] * v[2228] + v[453] * v[5220];
	v[2041] = v[161] * v[1795] + v[160] * v[1796] + v[451] * v[4915];
	v[2042] = v[158] * v[1795] + v[157] * v[1796] + v[445] * v[5211];
	v[2049] = v[155] * v[1821] + v[154] * v[1822] + v[427] * v[5208];
	v[2050] = v[152] * v[1821] + v[151] * v[1822] + v[1938] * v[4762] + v[425] * v[4916];
	v[2051] = v[149] * v[1821] + v[148] * v[1822] + v[419] * v[5200];
	v[2058] = v[155] * v[1809] + v[154] * v[1810] + v[1937] * v[4760] + v[427] * v[4917];
	v[2059] = v[152] * v[1809] + v[151] * v[1810] + v[425] * v[5204];
	v[2060] = v[149] * v[1809] + v[148] * v[1810] + v[419] * v[5199];
	v[2067] = v[155] * v[1797] + v[154] * v[1798] + v[1926] * v[2222] + v[427] * v[5207];
	v[2068] = v[152] * v[1797] + v[151] * v[1798] + v[425] * v[4918];
	v[2069] = v[149] * v[1797] + v[148] * v[1798] + v[419] * v[5198];
	v[2070] = -(v[1903] * v[281]) + v[1918] * v[297] - v[1914] * v[298] + v[1909] * v[4993];
	v[2071] = -(v[1919] * v[281]) - v[1904] * v[297] - v[1914] * v[299] + v[1902] * v[4992];
	v[2072] = -(v[1915] * v[281]) - v[1904] * v[298] + v[1918] * v[299] + v[1895] * v[4991];
	v[8630] = 0e0;
	v[8631] = 0e0;
	v[8632] = 0e0;
	v[8633] = 0e0;
	v[8634] = 0e0;
	v[8635] = 0e0;
	v[8636] = 0e0;
	v[8637] = 0e0;
	v[8638] = 0e0;
	v[8639] = 0e0;
	v[8640] = 0e0;
	v[8641] = 0e0;
	v[8642] = 0e0;
	v[8643] = 0e0;
	v[8644] = 0e0;
	v[8645] = 0e0;
	v[8646] = 0e0;
	v[8647] = 0e0;
	v[8648] = 0e0;
	v[8649] = 0e0;
	v[8650] = 0e0;
	v[8651] = v[2070];
	v[8652] = -v[2071];
	v[8653] = v[2072];
	v[2074] = v[2023] * v[281];
	v[2075] = v[2024] * v[281];
	v[2076] = v[2031] * v[281];
	v[2077] = v[2032] * v[4741];
	v[2078] = v[2033] * v[281];
	v[2079] = v[2040] * v[281];
	v[2080] = v[2041] * v[281];
	v[2081] = v[2042] * v[4741];
	v[2082] = 1e0 / (v[442] * v[442]);
	v[4928] = -(v[2082] * v[451]);
	v[4927] = -(v[2082] * v[453]);
	v[4926] = -(v[2082] * v[445]);
	v[4925] = -(v[2082] * v[440]);
	v[4924] = -(v[1737] * v[2082]);
	v[4923] = -(v[2082] * v[441]);
	v[4922] = -(v[2082] * v[4764]);
	v[3894] = -(v[2082] * v[5119]);
	v[3893] = -(v[2082] * v[5118]);
	v[3892] = -(v[2082] * v[5116]);
	v[3891] = -(v[2082] * v[5117]);
	v[3890] = -(v[2082] * v[5114]);
	v[3889] = -(v[2082] * v[5115]);
	v[3888] = -(v[2082] * v[443]);
	v[5122] = v[3888] * v[445];
	v[3887] = -(v[2082] * v[446]);
	v[5121] = v[3887] * v[451];
	v[5120] = -(v[2082] * v[453] * v[455]);
	v[3885] = -(v[2082] * v[4919]);
	v[3884] = -(v[2082] * v[4920]);
	v[3883] = -(v[2082] * v[4921]);
	v[3882] = v[1911] * v[4925];
	v[3881] = v[1910] * v[4923];
	v[3880] = v[1899] * v[4922];
	v[3664] = v[1738] * v[4922];
	v[3663] = v[4766] * v[4924];
	v[3661] = v[1739] * v[3888];
	v[5422] = v[3880] + (v[3661] + v[3663]) * v[453];
	v[5095] = v[3661] + v[3664];
	v[5423] = v[3663] + v[5095];
	v[3654] = v[1739] * v[4923];
	v[3651] = v[1738] * v[3887];
	v[5421] = v[3651] + v[3654];
	v[3650] = v[4768] * v[4924];
	v[5094] = v[3650] + v[3651];
	v[5420] = v[3654] + v[5094];
	v[5419] = v[3881] + v[453] * v[5094];
	v[3643] = v[1739] * v[4925];
	v[3641] = -(v[1738] * v[2082] * v[452]);
	v[3640] = v[455] * v[4924];
	v[5418] = v[3640] + v[3643];
	v[5093] = v[3640] + v[3641];
	v[5417] = v[3882] + v[451] * v[5093];
	v[5416] = v[3643] + v[5093];
	v[3506] = v[293] * v[4926];
	v[3502] = v[288] * v[4926];
	v[3499] = v[286] * v[4927];
	v[3498] = v[285] * v[4928];
	v[3493] = v[291] * v[4927];
	v[3492] = v[290] * v[4928];
	v[5396] = v[3492] + v[3502];
	v[5054] = v[3492] + v[3493];
	v[5394] = v[3502] + v[5054];
	v[3490] = v[284] * v[4926];
	v[5392] = v[3490] + v[3499];
	v[5055] = v[3490] + v[3498];
	v[5395] = v[3499] + v[5055];
	v[3485] = v[296] * v[4927];
	v[5397] = v[3485] + v[3506];
	v[3484] = v[295] * v[4928];
	v[5052] = v[3484] + v[3485];
	v[5393] = v[3506] + v[5052];
	v[2818] = v[3890] * v[445] + v[3889] * v[451] + v[5120];
	v[2813] = v[3892] * v[445] + v[3891] * v[453] + v[5121];
	v[2808] = v[3893] * v[451] + v[3894] * v[453] + v[5122];
	v[2616] = v[286] * v[4922];
	v[2614] = v[291] * v[4923];
	v[2610] = v[295] * v[4925];
	v[3879] = v[1895] + v[1902] + v[1909] + v[1911] * v[2610] + v[1910] * v[2614] + v[1899] * v[2616] + v[1739] * v[2808]
		+ v[1738] * v[2813] + v[1737] * v[2818];
	v[3877] = -(v[1915] * v[297]) - v[1919] * v[298] - v[1903] * v[299] + v[2022] * v[4784] + v[2032] * v[4785]
		+ v[2042] * v[4786] + v[3879] * v[5113] + v[2023] * v[722] + v[2024] * v[723] + v[2031] * v[724] + v[2033] * v[731]
		+ v[2040] * v[732] + v[2041] * v[737] - v[2070] * v[964] - v[2072] * v[965] + v[2071] * v[966];
	v[3997] = -v[2077] - v[2081] + v[3877] * v[749];
	v[3995] = v[2077] + v[3997] - v[2022] * v[4741];
	v[3993] = -v[2077] + v[2081] + v[3995];
	v[2084] = -(v[1930] * v[262]) + v[1945] * v[278] - v[1941] * v[279] + v[1936] * v[4989];
	v[2085] = -(v[1946] * v[262]) - v[1931] * v[278] - v[1941] * v[280] + v[1929] * v[4988];
	v[2086] = -(v[1942] * v[262]) - v[1931] * v[279] + v[1945] * v[280] + v[1922] * v[4987];
	v[8606] = 0e0;
	v[8607] = 0e0;
	v[8608] = 0e0;
	v[8609] = 0e0;
	v[8610] = 0e0;
	v[8611] = 0e0;
	v[8612] = 0e0;
	v[8613] = 0e0;
	v[8614] = 0e0;
	v[8615] = 0e0;
	v[8616] = 0e0;
	v[8617] = 0e0;
	v[8618] = 0e0;
	v[8619] = 0e0;
	v[8620] = 0e0;
	v[8621] = v[2084];
	v[8622] = -v[2085];
	v[8623] = v[2086];
	v[8624] = 0e0;
	v[8625] = 0e0;
	v[8626] = 0e0;
	v[8627] = 0e0;
	v[8628] = 0e0;
	v[8629] = 0e0;
	v[2088] = v[2050] * v[262];
	v[2089] = v[2051] * v[262];
	v[2090] = v[2058] * v[262];
	v[2091] = v[2059] * v[4737];
	v[2092] = v[2060] * v[262];
	v[2093] = v[2067] * v[262];
	v[2094] = v[2068] * v[262];
	v[2095] = v[2069] * v[4737];
	v[2096] = 1e0 / (v[416] * v[416]);
	v[4938] = -(v[2096] * v[425]);
	v[4937] = -(v[2096] * v[427]);
	v[4936] = -(v[2096] * v[419]);
	v[4935] = -(v[2096] * v[414]);
	v[4934] = -(v[1740] * v[2096]);
	v[4933] = -(v[2096] * v[415]);
	v[4932] = -(v[2096] * v[4759]);
	v[3916] = -(v[2096] * v[5132]);
	v[3915] = -(v[2096] * v[5131]);
	v[3914] = -(v[2096] * v[5129]);
	v[3913] = -(v[2096] * v[5130]);
	v[3912] = -(v[2096] * v[5127]);
	v[3911] = -(v[2096] * v[5128]);
	v[3910] = -(v[2096] * v[417]);
	v[5135] = v[3910] * v[419];
	v[3909] = -(v[2096] * v[420]);
	v[5134] = v[3909] * v[425];
	v[5133] = -(v[2096] * v[427] * v[429]);
	v[3907] = -(v[2096] * v[4929]);
	v[3906] = -(v[2096] * v[4930]);
	v[3905] = -(v[2096] * v[4931]);
	v[3904] = v[1938] * v[4935];
	v[3903] = v[1937] * v[4933];
	v[3902] = v[1926] * v[4932];
	v[3694] = v[1741] * v[4932];
	v[3693] = v[4761] * v[4934];
	v[3691] = v[1742] * v[3910];
	v[5430] = v[3902] + (v[3691] + v[3693]) * v[427];
	v[5098] = v[3691] + v[3694];
	v[5431] = v[3693] + v[5098];
	v[3684] = v[1742] * v[4933];
	v[3681] = v[1741] * v[3909];
	v[5429] = v[3681] + v[3684];
	v[3680] = v[4763] * v[4934];
	v[5097] = v[3680] + v[3681];
	v[5428] = v[3684] + v[5097];
	v[5427] = v[3903] + v[427] * v[5097];
	v[3673] = v[1742] * v[4935];
	v[3671] = -(v[1741] * v[2096] * v[426]);
	v[3670] = v[429] * v[4934];
	v[5426] = v[3670] + v[3673];
	v[5096] = v[3670] + v[3671];
	v[5425] = v[3904] + v[425] * v[5096];
	v[5424] = v[3673] + v[5096];
	v[3539] = v[274] * v[4936];
	v[3535] = v[269] * v[4936];
	v[3532] = v[267] * v[4937];
	v[3531] = v[266] * v[4938];
	v[3526] = v[272] * v[4937];
	v[3525] = v[271] * v[4938];
	v[5402] = v[3525] + v[3535];
	v[5062] = v[3525] + v[3526];
	v[5400] = v[3535] + v[5062];
	v[3523] = v[265] * v[4936];
	v[5398] = v[3523] + v[3532];
	v[5063] = v[3523] + v[3531];
	v[5401] = v[3532] + v[5063];
	v[3518] = v[277] * v[4937];
	v[5403] = v[3518] + v[3539];
	v[3517] = v[276] * v[4938];
	v[5060] = v[3517] + v[3518];
	v[5399] = v[3539] + v[5060];
	v[2803] = v[3912] * v[419] + v[3911] * v[425] + v[5133];
	v[2798] = v[3914] * v[419] + v[3913] * v[427] + v[5134];
	v[2793] = v[3915] * v[425] + v[3916] * v[427] + v[5135];
	v[2588] = v[267] * v[4932];
	v[2586] = v[272] * v[4933];
	v[2582] = v[276] * v[4935];
	v[3901] = v[1922] + v[1929] + v[1936] + v[1938] * v[2582] + v[1937] * v[2586] + v[1926] * v[2588] + v[1742] * v[2793]
		+ v[1741] * v[2798] + v[1740] * v[2803];
	v[3899] = -(v[1942] * v[278]) - v[1946] * v[279] - v[1930] * v[280] + v[2049] * v[4787] + v[2059] * v[4788]
		+ v[2069] * v[4789] + v[3901] * v[5126] + v[2050] * v[753] + v[2051] * v[754] + v[2058] * v[755] + v[2060] * v[762]
		+ v[2067] * v[763] + v[2068] * v[768] - v[2084] * v[970] - v[2086] * v[971] + v[2085] * v[972];
	v[3991] = -v[2091] - v[2095] + v[3899] * v[780];
	v[3989] = v[2091] + v[3991] - v[2049] * v[4737];
	v[3987] = -v[2091] + v[2095] + v[3989];
	v[2098] = v[2099] + v[1892] * v[239];
	v[2100] = -v[2099] + v[1892] * v[240];
	v[1800] = v[1800] + v[2098];
	v[1799] = v[1799] + v[2100];
	v[2101] = v[2102] + v[1891] * v[239];
	v[2103] = -v[2102] + v[1891] * v[240];
	v[1812] = v[1812] + v[2101];
	v[1811] = v[1811] + v[2103];
	v[2104] = v[2105] + v[1889] * v[239];
	v[2106] = -v[2105] + v[1889] * v[240];
	v[1824] = v[1824] + v[2104];
	v[1823] = v[1823] + v[2106];
	v[1825] = v[1825] + v[2106] * v[252] + v[4939] * v[509];
	v[1826] = v[1826] + v[2106] * v[250] + v[4939] * v[508];
	v[1827] = v[1827] + v[2104] * v[252] + v[4940] * v[509];
	v[1828] = v[1828] + v[2104] * v[250] + v[4940] * v[508];
	v[1813] = v[1813] + v[2103] * v[252] + v[4941] * v[509];
	v[1814] = v[1814] + v[2103] * v[250] + v[4941] * v[508];
	v[1815] = v[1815] + v[2101] * v[252] + v[4942] * v[509];
	v[1816] = v[1816] + v[2101] * v[250] + v[4942] * v[508];
	v[1801] = v[1801] + v[2100] * v[252] + v[4943] * v[509];
	v[1802] = v[1802] + v[2100] * v[250] + v[4943] * v[508];
	v[1803] = v[1803] + v[2098] * v[252] + v[4944] * v[509];
	v[1804] = v[1804] + v[2098] * v[250] + v[4944] * v[508];
	v[2119] = v[101] * v[1825] + v[100] * v[1826] + v[401] * v[5195];
	v[2120] = v[1965] * v[4757] + v[399] * v[4945] + v[1826] * v[97] + v[1825] * v[98];
	v[2121] = v[393] * v[5187] + v[1826] * v[94] + v[1825] * v[95];
	v[2128] = v[101] * v[1813] + v[100] * v[1814] + v[1964] * v[4755] + v[401] * v[4946];
	v[2129] = v[399] * v[5191] + v[1814] * v[97] + v[1813] * v[98];
	v[2130] = v[393] * v[5186] + v[1814] * v[94] + v[1813] * v[95];
	v[2137] = v[101] * v[1801] + v[100] * v[1802] + v[1953] * v[2216] + v[401] * v[5194];
	v[2138] = v[399] * v[4947] + v[1802] * v[97] + v[1801] * v[98];
	v[2139] = v[393] * v[5185] + v[1802] * v[94] + v[1801] * v[95];
	v[2146] = v[375] * v[5182] + v[1828] * v[91] + v[1827] * v[92];
	v[2147] = v[1992] * v[4752] + v[373] * v[4948] + v[1828] * v[88] + v[1827] * v[89];
	v[2148] = v[367] * v[5174] + v[1828] * v[85] + v[1827] * v[86];
	v[2155] = v[1991] * v[4750] + v[375] * v[4949] + v[1816] * v[91] + v[1815] * v[92];
	v[2156] = v[373] * v[5178] + v[1816] * v[88] + v[1815] * v[89];
	v[2157] = v[367] * v[5173] + v[1816] * v[85] + v[1815] * v[86];
	v[2164] = v[1980] * v[2210] + v[375] * v[5181] + v[1804] * v[91] + v[1803] * v[92];
	v[2165] = v[373] * v[4950] + v[1804] * v[88] + v[1803] * v[89];
	v[2166] = v[367] * v[5172] + v[1804] * v[85] + v[1803] * v[86];
	v[2167] = -(v[194] * v[1957]) + v[1972] * v[210] - v[1968] * v[211] + v[1963] * v[4979];
	v[2168] = -(v[194] * v[1973]) - v[1958] * v[210] - v[1968] * v[212] + v[1956] * v[4978];
	v[2169] = -(v[194] * v[1969]) - v[1958] * v[211] + v[1972] * v[212] + v[1949] * v[4977];
	v[8582] = 0e0;
	v[8583] = 0e0;
	v[8584] = 0e0;
	v[8585] = 0e0;
	v[8586] = 0e0;
	v[8587] = 0e0;
	v[8588] = 0e0;
	v[8589] = 0e0;
	v[8590] = 0e0;
	v[8591] = v[2167];
	v[8592] = -v[2168];
	v[8593] = v[2169];
	v[8594] = 0e0;
	v[8595] = 0e0;
	v[8596] = 0e0;
	v[8597] = 0e0;
	v[8598] = 0e0;
	v[8599] = 0e0;
	v[8600] = 0e0;
	v[8601] = 0e0;
	v[8602] = 0e0;
	v[8603] = 0e0;
	v[8604] = 0e0;
	v[8605] = 0e0;
	v[2171] = v[194] * v[2120];
	v[2172] = v[194] * v[2121];
	v[2173] = v[194] * v[2128];
	v[2174] = v[2129] * v[4728];
	v[2175] = v[194] * v[2130];
	v[2176] = v[194] * v[2137];
	v[2177] = v[194] * v[2138];
	v[2178] = v[2139] * v[4728];
	v[2179] = 1e0 / (v[390] * v[390]);
	v[4960] = -(v[2179] * v[399]);
	v[4959] = -(v[2179] * v[401]);
	v[4958] = -(v[2179] * v[393]);
	v[4957] = -(v[2179] * v[388]);
	v[4956] = -(v[1743] * v[2179]);
	v[4955] = -(v[2179] * v[389]);
	v[4954] = -(v[2179] * v[4754]);
	v[3938] = -(v[2179] * v[5145]);
	v[3937] = -(v[2179] * v[5144]);
	v[3936] = -(v[2179] * v[5142]);
	v[3935] = -(v[2179] * v[5143]);
	v[3934] = -(v[2179] * v[5140]);
	v[3933] = -(v[2179] * v[5141]);
	v[3932] = -(v[2179] * v[391]);
	v[5148] = v[393] * v[3932];
	v[3931] = -(v[2179] * v[394]);
	v[5147] = v[3931] * v[399];
	v[5146] = -(v[2179] * v[401] * v[403]);
	v[3929] = -(v[2179] * v[4951]);
	v[3928] = -(v[2179] * v[4952]);
	v[3927] = -(v[2179] * v[4953]);
	v[3926] = v[1965] * v[4957];
	v[3925] = v[1964] * v[4955];
	v[3924] = v[1953] * v[4954];
	v[3798] = v[1744] * v[4954];
	v[3797] = v[4756] * v[4956];
	v[3795] = v[1745] * v[3932];
	v[5438] = v[3924] + (v[3795] + v[3797]) * v[401];
	v[5109] = v[3795] + v[3798];
	v[5439] = v[3797] + v[5109];
	v[3788] = v[1745] * v[4955];
	v[3785] = v[1744] * v[3931];
	v[5437] = v[3785] + v[3788];
	v[3784] = v[4758] * v[4956];
	v[5108] = v[3784] + v[3785];
	v[5436] = v[3788] + v[5108];
	v[5435] = v[3925] + v[401] * v[5108];
	v[3777] = v[1745] * v[4957];
	v[3775] = -(v[1744] * v[2179] * v[400]);
	v[3774] = v[403] * v[4956];
	v[5434] = v[3774] + v[3777];
	v[5107] = v[3774] + v[3775];
	v[5433] = v[3926] + v[399] * v[5107];
	v[5432] = v[3777] + v[5107];
	v[3572] = v[206] * v[4958];
	v[3568] = v[201] * v[4958];
	v[3565] = v[199] * v[4959];
	v[3564] = v[198] * v[4960];
	v[3559] = v[204] * v[4959];
	v[3558] = v[203] * v[4960];
	v[5408] = v[3558] + v[3568];
	v[5070] = v[3558] + v[3559];
	v[5406] = v[3568] + v[5070];
	v[3556] = v[197] * v[4958];
	v[5404] = v[3556] + v[3565];
	v[5071] = v[3556] + v[3564];
	v[5407] = v[3565] + v[5071];
	v[3551] = v[209] * v[4959];
	v[5409] = v[3551] + v[3572];
	v[3550] = v[208] * v[4960];
	v[5068] = v[3550] + v[3551];
	v[5405] = v[3572] + v[5068];
	v[2788] = v[393] * v[3934] + v[3933] * v[399] + v[5146];
	v[2783] = v[393] * v[3936] + v[3935] * v[401] + v[5147];
	v[2778] = v[3937] * v[399] + v[3938] * v[401] + v[5148];
	v[2560] = v[199] * v[4954];
	v[2558] = v[204] * v[4955];
	v[2554] = v[208] * v[4957];
	v[3923] = v[1949] + v[1956] + v[1963] + v[1965] * v[2554] + v[1964] * v[2558] + v[1953] * v[2560] + v[1745] * v[2778]
		+ v[1744] * v[2783] + v[1743] * v[2788];
	v[3921] = -(v[1969] * v[210]) - v[1973] * v[211] - v[1957] * v[212] + v[2119] * v[4796] + v[2129] * v[4797]
		+ v[2139] * v[4798] + v[3923] * v[5139] + v[2120] * v[904] + v[2121] * v[905] + v[2128] * v[906] + v[2130] * v[913]
		+ v[2137] * v[914] + v[2138] * v[919] - v[2167] * v[976] - v[2169] * v[977] + v[2168] * v[978];
	v[3985] = -v[2174] - v[2178] + v[3921] * v[931];
	v[3983] = v[2174] + v[3985] - v[2119] * v[4728];
	v[3981] = -v[2174] + v[2178] + v[3983];
	v[2181] = -(v[175] * v[1984]) - v[192] * v[1995] + v[191] * v[1999] + v[1990] * v[4975];
	v[2182] = -(v[191] * v[1985]) - v[193] * v[1995] - v[175] * v[2000] + v[1983] * v[4974];
	v[2183] = -(v[192] * v[1985]) - v[175] * v[1996] + v[193] * v[1999] + v[1976] * v[4973];
	v[8558] = 0e0;
	v[8559] = 0e0;
	v[8560] = 0e0;
	v[8561] = v[2181];
	v[8562] = -v[2182];
	v[8563] = v[2183];
	v[8564] = 0e0;
	v[8565] = 0e0;
	v[8566] = 0e0;
	v[8567] = 0e0;
	v[8568] = 0e0;
	v[8569] = 0e0;
	v[8570] = 0e0;
	v[8571] = 0e0;
	v[8572] = 0e0;
	v[8573] = 0e0;
	v[8574] = 0e0;
	v[8575] = 0e0;
	v[8576] = 0e0;
	v[8577] = 0e0;
	v[8578] = 0e0;
	v[8579] = 0e0;
	v[8580] = 0e0;
	v[8581] = 0e0;
	v[2185] = v[175] * v[2147];
	v[2186] = v[175] * v[2148];
	v[2187] = v[175] * v[2155];
	v[2188] = v[2156] * v[4724];
	v[2189] = v[175] * v[2157];
	v[2190] = v[175] * v[2164];
	v[2191] = v[175] * v[2165];
	v[2192] = v[2166] * v[4724];
	v[2193] = 1e0 / (v[364] * v[364]);
	v[4970] = -(v[2193] * v[373]);
	v[4969] = -(v[2193] * v[375]);
	v[4968] = -(v[2193] * v[367]);
	v[4967] = -(v[2193] * v[362]);
	v[4966] = -(v[1746] * v[2193]);
	v[4965] = -(v[2193] * v[363]);
	v[4964] = -(v[2193] * v[4749]);
	v[3960] = -(v[2193] * v[5158]);
	v[3959] = -(v[2193] * v[5157]);
	v[3958] = -(v[2193] * v[5155]);
	v[3957] = -(v[2193] * v[5156]);
	v[3956] = -(v[2193] * v[5153]);
	v[3955] = -(v[2193] * v[5154]);
	v[3954] = -(v[2193] * v[365]);
	v[5161] = v[367] * v[3954];
	v[3953] = -(v[2193] * v[368]);
	v[5160] = v[373] * v[3953];
	v[5159] = -(v[2193] * v[375] * v[377]);
	v[3951] = -(v[2193] * v[4961]);
	v[3950] = -(v[2193] * v[4962]);
	v[3949] = -(v[2193] * v[4963]);
	v[3948] = v[1992] * v[4967];
	v[3947] = v[1991] * v[4965];
	v[3946] = v[1980] * v[4964];
	v[3828] = v[1747] * v[4964];
	v[3827] = v[4751] * v[4966];
	v[3825] = v[1748] * v[3954];
	v[5446] = v[375] * (v[3825] + v[3827]) + v[3946];
	v[5112] = v[3825] + v[3828];
	v[5447] = v[3827] + v[5112];
	v[3818] = v[1748] * v[4965];
	v[3815] = v[1747] * v[3953];
	v[5445] = v[3815] + v[3818];
	v[3814] = v[4753] * v[4966];
	v[5111] = v[3814] + v[3815];
	v[5444] = v[3818] + v[5111];
	v[5443] = v[3947] + v[375] * v[5111];
	v[3807] = v[1748] * v[4967];
	v[3805] = -(v[1747] * v[2193] * v[374]);
	v[3804] = v[377] * v[4966];
	v[5442] = v[3804] + v[3807];
	v[5110] = v[3804] + v[3805];
	v[5441] = v[3948] + v[373] * v[5110];
	v[5440] = v[3807] + v[5110];
	v[3605] = v[187] * v[4968];
	v[3601] = v[182] * v[4968];
	v[3598] = v[180] * v[4969];
	v[3597] = v[179] * v[4970];
	v[3592] = v[185] * v[4969];
	v[3591] = v[184] * v[4970];
	v[5414] = v[3591] + v[3601];
	v[5078] = v[3591] + v[3592];
	v[5412] = v[3601] + v[5078];
	v[3589] = v[178] * v[4968];
	v[5410] = v[3589] + v[3598];
	v[5079] = v[3589] + v[3597];
	v[5413] = v[3598] + v[5079];
	v[3584] = v[190] * v[4969];
	v[5415] = v[3584] + v[3605];
	v[3583] = v[189] * v[4970];
	v[5076] = v[3583] + v[3584];
	v[5411] = v[3605] + v[5076];
	v[2773] = v[373] * v[3955] + v[367] * v[3956] + v[5159];
	v[2768] = v[375] * v[3957] + v[367] * v[3958] + v[5160];
	v[2763] = v[373] * v[3959] + v[375] * v[3960] + v[5161];
	v[2532] = v[180] * v[4964];
	v[2530] = v[185] * v[4965];
	v[2526] = v[189] * v[4967];
	v[3945] = v[1976] + v[1983] + v[1990] + v[1992] * v[2526] + v[1991] * v[2530] + v[1980] * v[2532] + v[1748] * v[2763]
		+ v[1747] * v[2768] + v[1746] * v[2773];
	v[3943] = -(v[193] * v[1984]) - v[191] * v[1996] - v[192] * v[2000] + v[2146] * v[4799] + v[2156] * v[4800]
		+ v[2166] * v[4801] + v[3945] * v[5152] + v[2147] * v[935] + v[2148] * v[936] + v[2155] * v[937] + v[2157] * v[944]
		+ v[2164] * v[945] + v[2165] * v[950] - v[2181] * v[982] - v[2183] * v[983] + v[2182] * v[984];
	v[3979] = -v[2188] - v[2192] + v[3943] * v[962];
	v[3977] = v[2188] + v[3979] - v[2146] * v[4724];
	v[3975] = -v[2188] + v[2192] + v[3977];
	v[2195] = v[2075] + v[2079];
	v[2196] = v[2074] + v[2076];
	v[2197] = v[2078] + v[2080];
	v[2198] = v[2089] + v[2093];
	v[2199] = v[2088] + v[2090];
	v[2200] = v[2092] + v[2094];
	v[2201] = v[2172] + v[2176];
	v[2202] = v[2171] + v[2173];
	v[2203] = v[2175] + v[2177];
	v[2204] = v[2186] + v[2190];
	v[2205] = v[2185] + v[2187];
	v[2206] = v[2189] + v[2191];
	v[6127] = -(v[1498] * v[1566]) - v[1499] * v[1567] - v[1500] * v[1568] + v[1800] + (-(v[1566] * v[1665])
		- v[1567] * v[1666] - v[1568] * v[1667]) * v[4971];
	v[6128] = -(v[1498] * v[1570]) - v[1499] * v[1571] - v[1500] * v[1572] + v[1812] + (-(v[1570] * v[1665])
		- v[1571] * v[1666] - v[1572] * v[1667]) * v[4971];
	v[6129] = -(v[1498] * v[1574]) - v[1499] * v[1575] - v[1500] * v[1576] + v[1824] + (-(v[1574] * v[1665])
		- v[1575] * v[1666] - v[1576] * v[1667]) * v[4971];
	v[6130] = -(v[1498] * v[1578]) - v[1499] * v[1579] - v[1500] * v[1580] + v[2185] - v[2187] - v[2181] * v[4724] + v[19] *
		(v[2159] + v[1747] * v[4750] + v[1746] * v[4752]) + v[2204] * v[983] + v[2206] * v[984] + v[3975] * v[986];
	v[6131] = -(v[1498] * v[1582]) - v[1499] * v[1583] - v[1500] * v[1584] - v[2186] + v[2190] + v[19] * (v[2151]
		+ v[1746] * v[2209] + v[1748] * v[2210]) + v[2182] * v[4724] + v[2206] * v[982] + v[2205] * v[983] + v[3977] * v[985];
	v[6132] = -(v[1498] * v[1586]) - v[1499] * v[1587] - v[1500] * v[1588] + v[2189] - v[2191] + v[19] * (v[2141]
		+ v[1747] * v[2211] + v[1748] * v[2212]) - v[2183] * v[4724] + v[3979] * v[981] + v[2204] * v[982] + v[2205] * v[984];
	v[6133] = -(v[1498] * v[1590]) - v[1499] * v[1591] - v[1500] * v[1592] + v[1799] + (-(v[1590] * v[1665])
		- v[1591] * v[1666] - v[1592] * v[1667]) * v[4971];
	v[6134] = -(v[1498] * v[1594]) - v[1499] * v[1595] - v[1500] * v[1596] + v[1811] + (-(v[1594] * v[1665])
		- v[1595] * v[1666] - v[1596] * v[1667]) * v[4971];
	v[6135] = -(v[1498] * v[1598]) - v[1499] * v[1599] - v[1500] * v[1600] + v[1823] + (-(v[1598] * v[1665])
		- v[1599] * v[1666] - v[1600] * v[1667]) * v[4971];
	v[6136] = -(v[1498] * v[1602]) - v[1499] * v[1603] - v[1500] * v[1604] + v[2171] - v[2173] - v[2167] * v[4728] + v[19] *
		(v[2132] + v[1744] * v[4755] + v[1743] * v[4757]) + v[2201] * v[977] + v[2203] * v[978] + v[3981] * v[980];
	v[6137] = -(v[1498] * v[1606]) - v[1499] * v[1607] - v[1500] * v[1608] - v[2172] + v[2176] + v[19] * (v[2124]
		+ v[1743] * v[2215] + v[1745] * v[2216]) + v[2168] * v[4728] + v[2203] * v[976] + v[2202] * v[977] + v[3983] * v[979];
	v[6138] = -(v[1498] * v[1610]) - v[1499] * v[1611] - v[1500] * v[1612] + v[2175] - v[2177] + v[19] * (v[2114]
		+ v[1744] * v[2217] + v[1745] * v[2218]) - v[2169] * v[4728] + v[3985] * v[975] + v[2201] * v[976] + v[2202] * v[978];
	v[6139] = -(v[1498] * v[1614]) - v[1499] * v[1615] - v[1500] * v[1616] + v[1794] + (-(v[1614] * v[1665])
		- v[1615] * v[1666] - v[1616] * v[1667]) * v[4971];
	v[6140] = -(v[1498] * v[1618]) - v[1499] * v[1619] - v[1500] * v[1620] + v[1806] + (-(v[1618] * v[1665])
		- v[1619] * v[1666] - v[1620] * v[1667]) * v[4971];
	v[6141] = -(v[1498] * v[1622]) - v[1499] * v[1623] - v[1500] * v[1624] + v[1818] + (-(v[1622] * v[1665])
		- v[1623] * v[1666] - v[1624] * v[1667]) * v[4971];
	v[6142] = -(v[1498] * v[1626]) - v[1499] * v[1627] - v[1500] * v[1628] + v[2088] - v[2090] - v[2084] * v[4737] + v[19] *
		(v[2062] + v[1741] * v[4760] + v[1740] * v[4762]) + v[2198] * v[971] + v[2200] * v[972] + v[3987] * v[974];
	v[6143] = -(v[1498] * v[1630]) - v[1499] * v[1631] - v[1500] * v[1632] - v[2089] + v[2093] + v[19] * (v[2054]
		+ v[1740] * v[2221] + v[1742] * v[2222]) + v[2085] * v[4737] + v[2200] * v[970] + v[2199] * v[971] + v[3989] * v[973];
	v[6144] = -(v[1498] * v[1634]) - v[1499] * v[1635] - v[1500] * v[1636] + v[2092] - v[2094] + v[19] * (v[2044]
		+ v[1741] * v[2223] + v[1742] * v[2224]) - v[2086] * v[4737] + v[3991] * v[969] + v[2198] * v[970] + v[2199] * v[972];
	v[6145] = -(v[1498] * v[1638]) - v[1499] * v[1639] - v[1500] * v[1640] + v[1793] + (-(v[1638] * v[1665])
		- v[1639] * v[1666] - v[1640] * v[1667]) * v[4971];
	v[6146] = -(v[1498] * v[1642]) - v[1499] * v[1643] - v[1500] * v[1644] + v[1805] + (-(v[1642] * v[1665])
		- v[1643] * v[1666] - v[1644] * v[1667]) * v[4971];
	v[6147] = -(v[1498] * v[1646]) - v[1499] * v[1647] - v[1500] * v[1648] + v[1817] + (-(v[1646] * v[1665])
		- v[1647] * v[1666] - v[1648] * v[1667]) * v[4971];
	v[6148] = -(v[1498] * v[1650]) - v[1499] * v[1651] - v[1500] * v[1652] + v[2074] - v[2076] - v[2070] * v[4741] + v[19] *
		(v[2035] + v[1738] * v[4765] + v[1737] * v[4767]) + v[2195] * v[965] + v[2197] * v[966] + v[3993] * v[968];
	v[6149] = -(v[1498] * v[1654]) - v[1499] * v[1655] - v[1500] * v[1656] - v[2075] + v[2079] + v[19] * (v[2027]
		+ v[1737] * v[2227] + v[1739] * v[2228]) + v[2071] * v[4741] + v[2197] * v[964] + v[2196] * v[965] + v[3995] * v[967];
	v[6150] = -(v[1498] * v[1658]) - v[1499] * v[1659] - v[1500] * v[1660] + v[2078] - v[2080] + v[19] * (v[2017]
		+ v[1738] * v[2229] + v[1739] * v[2230]) - v[2072] * v[4741] + v[3997] * v[963] + v[2195] * v[964] + v[2196] * v[966];
	for (i1663 = 1; i1663 <= 24; i1663++) {
		i4999 = (i1663 == 13 ? 1 : 0);
		i4998 = (i1663 == 19 ? 1 : 0);
		i4997 = (i1663 == 14 ? 1 : 0);
		i4996 = (i1663 == 20 ? 1 : 0);
		i4995 = (i1663 == 15 ? 1 : 0);
		i4994 = (i1663 == 21 ? 1 : 0);
		i4985 = (i1663 == 1 ? 1 : 0);
		i4984 = (i1663 == 7 ? 1 : 0);
		i4983 = (i1663 == 2 ? 1 : 0);
		i4982 = (i1663 == 8 ? 1 : 0);
		i4981 = (i1663 == 3 ? 1 : 0);
		i4980 = (i1663 == 9 ? 1 : 0);
		v[2237] = v[6157 + i1663];
		v[2238] = v[6181 + i1663];
		v[2239] = v[6205 + i1663];
		v[2240] = v[6229 + i1663];
		v[2241] = v[6253 + i1663];
		v[2242] = v[6277 + i1663];
		v[2243] = v[6301 + i1663];
		v[2244] = v[6325 + i1663];
		v[2245] = v[6349 + i1663];
		v[2246] = v[6373 + i1663];
		v[2247] = v[6397 + i1663];
		v[2248] = v[6421 + i1663];
		v[2249] = v[6445 + i1663];
		v[3944] = 2e0 * v[2249] * v[962];
		v[4972] = v[3944] / 2e0;
		v[3867] = v[3944] / 4e0;
		v[2329] = v[175] * v[3944];
		v[5000] = v[2329] * v[364];
		v[2250] = v[6469 + i1663];
		v[2251] = v[6565 + i1663];
		v[2252] = v[6661 + i1663];
		v[2253] = v[6757 + i1663];
		v[3922] = 2e0 * v[2253] * v[931];
		v[4976] = v[3922] / 2e0;
		v[3845] = v[3922] / 4e0;
		v[2351] = v[194] * v[3922];
		v[5001] = v[2351] * v[390];
		v[2254] = v[6781 + i1663];
		v[2255] = v[6877 + i1663];
		v[2256] = v[6973 + i1663];
		v[2263] = v[7213 + i1663];
		v[3900] = 2e0 * v[2263] * v[780];
		v[4986] = v[3900] / 2e0;
		v[3733] = v[3900] / 4e0;
		v[2433] = v[262] * v[3900];
		v[5002] = v[2433] * v[416];
		v[2264] = v[7237 + i1663];
		v[2265] = v[7333 + i1663];
		v[2266] = v[7429 + i1663];
		v[2267] = v[7525 + i1663];
		v[3878] = 2e0 * v[2267] * v[749];
		v[4990] = v[3878] / 2e0;
		v[3711] = v[3878] / 4e0;
		v[2455] = v[281] * v[3878];
		v[5003] = v[2455] * v[442];
		v[2268] = v[7549 + i1663];
		v[2269] = v[7645 + i1663];
		v[2270] = v[7741 + i1663];
		v[2277] = (i1663 == 24 ? 1 : 0);
		v[5051] = v[19] * v[2277];
		v[2278] = (i1663 == 23 ? 1 : 0);
		v[5049] = v[19] * v[2278];
		v[2279] = (i1663 == 22 ? 1 : 0);
		v[5050] = v[19] * v[2279];
		v[2280] = (i1663 == 18 ? 1 : 0);
		v[5059] = v[19] * v[2280];
		v[2281] = (i1663 == 17 ? 1 : 0);
		v[5057] = v[19] * v[2281];
		v[2282] = (i1663 == 16 ? 1 : 0);
		v[5058] = v[19] * v[2282];
		v[2283] = (i1663 == 12 ? 1 : 0);
		v[5067] = v[19] * v[2283];
		v[2284] = (i1663 == 11 ? 1 : 0);
		v[5065] = v[19] * v[2284];
		v[2285] = (i1663 == 10 ? 1 : 0);
		v[5066] = v[19] * v[2285];
		v[2286] = (i1663 == 6 ? 1 : 0);
		v[5075] = v[19] * v[2286];
		v[2287] = (i1663 == 5 ? 1 : 0);
		v[5073] = v[19] * v[2287];
		v[2288] = (i1663 == 4 ? 1 : 0);
		v[5074] = v[19] * v[2288];
		v[2289] = v[2237] - v[2286];
		v[2290] = v[2237] + v[2286];
		v[2291] = v[2238] - v[2288];
		v[2292] = v[2238] + v[2288];
		v[2293] = v[2239] + v[2287];
		v[2294] = v[2239] - v[2287];
		v[2295] = v[2240] - v[2283];
		v[2296] = v[2240] + v[2283];
		v[2297] = v[2241] - v[2285];
		v[2298] = v[2241] + v[2285];
		v[2299] = v[2242] + v[2284];
		v[2300] = v[2242] - v[2284];
		v[2301] = v[2243] - v[2280];
		v[2302] = v[2243] + v[2280];
		v[2303] = v[2244] - v[2282];
		v[2304] = v[2244] + v[2282];
		v[2305] = v[2245] + v[2281];
		v[2306] = v[2245] - v[2281];
		v[2307] = v[2246] - v[2277];
		v[2308] = v[2246] + v[2277];
		v[2309] = v[2247] - v[2279];
		v[2310] = v[2247] + v[2279];
		v[2311] = v[2248] + v[2278];
		v[2312] = v[2248] - v[2278];
		v[2313] = -(v[2286] * v[4724]) - v[4972] * v[983];
		v[2314] = v[2287] * v[4724] + v[4972] * v[984];
		v[2315] = -(v[2288] * v[4724]) - v[4972] * v[982];
		v[2316] = v[2250] * v[4724] + v[4801] * v[4972];
		v[2359] = v[2316] * v[367];
		v[2317] = v[175] * v[2289] + v[4972] * v[950];
		v[2318] = v[175] * v[2293] + v[4972] * v[945];
		v[2360] = v[2318] * v[375];
		v[2319] = v[175] * v[2290] + v[4972] * v[944];
		v[2320] = v[2251] * v[4724] + v[4800] * v[4972];
		v[2364] = v[2320] * v[373];
		v[2321] = v[175] * v[2291] + v[4972] * v[937];
		v[2365] = v[2321] * v[375];
		v[2322] = v[175] * v[2294] + v[4972] * v[936];
		v[2370] = v[2322] * v[367];
		v[2323] = v[175] * v[2292] + v[4972] * v[935];
		v[2324] = v[2252] * v[4724] + v[4799] * v[4972];
		v[2368] = v[2324] * v[375];
		v[2325] = -(v[175] * v[2313]) - v[191] * v[4972];
		v[2326] = v[2329] + v[2313] * v[4973];
		v[2327] = -(v[175] * v[2314]) - v[192] * v[4972];
		v[2328] = -(v[192] * v[2313]) - v[191] * v[2314];
		v[2330] = v[2329] + v[2314] * v[4974];
		v[2548] = v[1747] * v[2330];
		v[2331] = v[193] * v[2313] + v[191] * v[2315];
		v[2332] = -(v[193] * v[2314]) - v[192] * v[2315];
		v[2333] = v[2329] + v[2315] * v[4975];
		v[2544] = v[1748] * v[2333];
		v[2334] = -(v[175] * v[2315]) - v[193] * v[4972];
		v[2335] = -(v[2283] * v[4728]) - v[4976] * v[977];
		v[2336] = v[2284] * v[4728] + v[4976] * v[978];
		v[2337] = -(v[2285] * v[4728]) - v[4976] * v[976];
		v[2338] = v[2254] * v[4728] + v[4798] * v[4976];
		v[2374] = v[2338] * v[393];
		v[2339] = v[194] * v[2295] + v[4976] * v[919];
		v[2340] = v[194] * v[2299] + v[4976] * v[914];
		v[2375] = v[2340] * v[401];
		v[2341] = v[194] * v[2296] + v[4976] * v[913];
		v[2342] = v[2255] * v[4728] + v[4797] * v[4976];
		v[2379] = v[2342] * v[399];
		v[2343] = v[194] * v[2297] + v[4976] * v[906];
		v[2380] = v[2343] * v[401];
		v[2344] = v[194] * v[2300] + v[4976] * v[905];
		v[2385] = v[2344] * v[393];
		v[2345] = v[194] * v[2298] + v[4976] * v[904];
		v[2346] = v[2256] * v[4728] + v[4796] * v[4976];
		v[2383] = v[2346] * v[401];
		v[2347] = -(v[194] * v[2335]) - v[210] * v[4976];
		v[2348] = v[2351] + v[2335] * v[4977];
		v[2349] = -(v[194] * v[2336]) - v[211] * v[4976];
		v[2350] = -(v[211] * v[2335]) - v[210] * v[2336];
		v[2352] = v[2351] + v[2336] * v[4978];
		v[2576] = v[1744] * v[2352];
		v[2353] = v[212] * v[2335] + v[210] * v[2337];
		v[2354] = -(v[212] * v[2336]) - v[211] * v[2337];
		v[2355] = v[2351] + v[2337] * v[4979];
		v[2572] = v[1745] * v[2355];
		v[2356] = -(v[194] * v[2337]) - v[212] * v[4976];
		v[2357] = v[2359] + v[2317] * v[373];
		v[2358] = v[2357] + v[2360] + v[5074];
		v[2361] = v[2359] + v[2360];
		v[2362] = v[2364] + v[2319] * v[367];
		v[2363] = v[2362] + v[2365] + v[5073];
		v[2366] = v[2364] + v[2365];
		v[2367] = v[2368] + v[2370];
		v[2369] = v[2368] + v[2323] * v[373];
		v[2371] = v[2369] + v[2370] + v[5075];
		v[2372] = v[2374] + v[2339] * v[399];
		v[2373] = v[2372] + v[2375] + v[5066];
		v[2376] = v[2374] + v[2375];
		v[2377] = v[2379] + v[2341] * v[393];
		v[2378] = v[2377] + v[2380] + v[5065];
		v[2381] = v[2379] + v[2380];
		v[2382] = v[2383] + v[2385];
		v[2384] = v[2383] + v[2345] * v[399];
		v[2386] = v[2384] + v[2385] + v[5067];
		v[2387] = i4985 * v[19];
		v[2388] = i4983 * v[19];
		v[2389] = i4981 * v[19];
		v[2390] = i4984 * v[19];
		v[2391] = i4982 * v[19];
		v[2392] = i4980 * v[19];
		v[2393] = v[2316] * v[85] + v[2317] * v[88] + v[2318] * v[91];
		v[2394] = v[2316] * v[86] + v[2317] * v[89] + v[2318] * v[92];
		v[2395] = v[100] * v[2340] + v[2338] * v[94] + v[2339] * v[97];
		v[2396] = v[101] * v[2340] + v[2338] * v[95] + v[2339] * v[98];
		v[2397] = v[2319] * v[85] + v[2320] * v[88] + v[2321] * v[91];
		v[2398] = v[2319] * v[86] + v[2320] * v[89] + v[2321] * v[92];
		v[2399] = v[100] * v[2343] + v[2341] * v[94] + v[2342] * v[97];
		v[2400] = v[101] * v[2343] + v[2341] * v[95] + v[2342] * v[98];
		v[2401] = v[2322] * v[85] + v[2323] * v[88] + v[2324] * v[91];
		v[2402] = v[2322] * v[86] + v[2323] * v[89] + v[2324] * v[92];
		v[2403] = v[100] * v[2346] + v[2344] * v[94] + v[2345] * v[97];
		v[2404] = v[101] * v[2346] + v[2344] * v[95] + v[2345] * v[98];
		v[2405] = i4980 + v[2403] * v[250] + v[2404] * v[252];
		v[2406] = i4980;
		v[2407] = i4981 + v[2401] * v[250] + v[2402] * v[252];
		v[2408] = i4981;
		v[2409] = i4982 + v[2399] * v[250] + v[2400] * v[252];
		v[2410] = i4982;
		v[2411] = i4983 + v[2397] * v[250] + v[2398] * v[252];
		v[2412] = i4983;
		v[2413] = i4984 + v[2395] * v[250] + v[2396] * v[252];
		v[2414] = i4984;
		v[2415] = i4985 + v[2393] * v[250] + v[2394] * v[252];
		v[2416] = i4985;
		v[2417] = -(v[2280] * v[4737]) - v[4986] * v[971];
		v[2418] = v[2281] * v[4737] + v[4986] * v[972];
		v[2419] = -(v[2282] * v[4737]) - v[4986] * v[970];
		v[2420] = v[2264] * v[4737] + v[4789] * v[4986];
		v[2463] = v[2420] * v[419];
		v[2421] = v[2301] * v[262] + v[4986] * v[768];
		v[2422] = v[2305] * v[262] + v[4986] * v[763];
		v[2464] = v[2422] * v[427];
		v[2423] = v[2302] * v[262] + v[4986] * v[762];
		v[2424] = v[2265] * v[4737] + v[4788] * v[4986];
		v[2468] = v[2424] * v[425];
		v[2425] = v[2303] * v[262] + v[4986] * v[755];
		v[2469] = v[2425] * v[427];
		v[2426] = v[2306] * v[262] + v[4986] * v[754];
		v[2474] = v[2426] * v[419];
		v[2427] = v[2304] * v[262] + v[4986] * v[753];
		v[2428] = v[2266] * v[4737] + v[4787] * v[4986];
		v[2472] = v[2428] * v[427];
		v[2429] = -(v[2417] * v[262]) - v[278] * v[4986];
		v[2430] = v[2433] + v[2417] * v[4987];
		v[2431] = -(v[2418] * v[262]) - v[279] * v[4986];
		v[2432] = -(v[2418] * v[278]) - v[2417] * v[279];
		v[2434] = v[2433] + v[2418] * v[4988];
		v[2604] = v[1741] * v[2434];
		v[2435] = v[2419] * v[278] + v[2417] * v[280];
		v[2436] = -(v[2419] * v[279]) - v[2418] * v[280];
		v[2437] = v[2433] + v[2419] * v[4989];
		v[2600] = v[1742] * v[2437];
		v[2438] = -(v[2419] * v[262]) - v[280] * v[4986];
		v[2439] = -(v[2277] * v[4741]) - v[4990] * v[965];
		v[2440] = v[2278] * v[4741] + v[4990] * v[966];
		v[2441] = -(v[2279] * v[4741]) - v[4990] * v[964];
		v[2442] = v[2268] * v[4741] + v[4786] * v[4990];
		v[2478] = v[2442] * v[445];
		v[2443] = v[2307] * v[281] + v[4990] * v[737];
		v[2444] = v[2311] * v[281] + v[4990] * v[732];
		v[2479] = v[2444] * v[453];
		v[2445] = v[2308] * v[281] + v[4990] * v[731];
		v[2446] = v[2269] * v[4741] + v[4785] * v[4990];
		v[2483] = v[2446] * v[451];
		v[2447] = v[2309] * v[281] + v[4990] * v[724];
		v[2484] = v[2447] * v[453];
		v[2448] = v[2312] * v[281] + v[4990] * v[723];
		v[2489] = v[2448] * v[445];
		v[2449] = v[2310] * v[281] + v[4990] * v[722];
		v[2450] = v[2270] * v[4741] + v[4784] * v[4990];
		v[2487] = v[2450] * v[453];
		v[2451] = -(v[2439] * v[281]) - v[297] * v[4990];
		v[2452] = v[2455] + v[2439] * v[4991];
		v[2453] = -(v[2440] * v[281]) - v[298] * v[4990];
		v[2454] = -(v[2440] * v[297]) - v[2439] * v[298];
		v[2456] = v[2455] + v[2440] * v[4992];
		v[2632] = v[1738] * v[2456];
		v[2457] = v[2441] * v[297] + v[2439] * v[299];
		v[2458] = -(v[2441] * v[298]) - v[2440] * v[299];
		v[2459] = v[2455] + v[2441] * v[4993];
		v[2628] = v[1739] * v[2459];
		v[2460] = -(v[2441] * v[281]) - v[299] * v[4990];
		v[2461] = v[2463] + v[2421] * v[425];
		v[2462] = v[2461] + v[2464] + v[5058];
		v[2465] = v[2463] + v[2464];
		v[2466] = v[2468] + v[2423] * v[419];
		v[2467] = v[2466] + v[2469] + v[5057];
		v[2470] = v[2468] + v[2469];
		v[2471] = v[2472] + v[2474];
		v[2473] = v[2472] + v[2427] * v[425];
		v[2475] = v[2473] + v[2474] + v[5059];
		v[2476] = v[2478] + v[2443] * v[451];
		v[2477] = v[2476] + v[2479] + v[5050];
		v[2480] = v[2478] + v[2479];
		v[2481] = v[2483] + v[2445] * v[445];
		v[2482] = v[2481] + v[2484] + v[5049];
		v[2485] = v[2483] + v[2484];
		v[2486] = v[2487] + v[2489];
		v[2488] = v[2487] + v[2449] * v[451];
		v[2490] = v[2488] + v[2489] + v[5051];
		v[2491] = i4999 * v[19];
		v[2492] = i4997 * v[19];
		v[2493] = i4995 * v[19];
		v[2494] = i4998 * v[19];
		v[2495] = i4996 * v[19];
		v[2496] = i4994 * v[19];
		v[2497] = v[148] * v[2420] + v[151] * v[2421] + v[154] * v[2422];
		v[2498] = v[149] * v[2420] + v[152] * v[2421] + v[155] * v[2422];
		v[2499] = v[157] * v[2442] + v[160] * v[2443] + v[163] * v[2444];
		v[2500] = v[158] * v[2442] + v[161] * v[2443] + v[164] * v[2444];
		v[2501] = v[148] * v[2423] + v[151] * v[2424] + v[154] * v[2425];
		v[2502] = v[149] * v[2423] + v[152] * v[2424] + v[155] * v[2425];
		v[2503] = v[157] * v[2445] + v[160] * v[2446] + v[163] * v[2447];
		v[2504] = v[158] * v[2445] + v[161] * v[2446] + v[164] * v[2447];
		v[2505] = v[148] * v[2426] + v[151] * v[2427] + v[154] * v[2428];
		v[2506] = v[149] * v[2426] + v[152] * v[2427] + v[155] * v[2428];
		v[2507] = v[157] * v[2448] + v[160] * v[2449] + v[163] * v[2450];
		v[2508] = v[158] * v[2448] + v[161] * v[2449] + v[164] * v[2450];
		v[2509] = i4994 + v[2507] * v[337] + v[2508] * v[339];
		v[2510] = i4994;
		v[2511] = i4995 + v[2505] * v[337] + v[2506] * v[339];
		v[2512] = i4995;
		v[2513] = v[240] * v[2405] + v[239] * v[2407] - v[2511] * v[326] - v[2509] * v[327];
		v[2514] = i4996 + v[2503] * v[337] + v[2504] * v[339];
		v[2515] = i4996;
		v[2516] = i4997 + v[2501] * v[337] + v[2502] * v[339];
		v[2517] = i4997;
		v[2518] = v[240] * v[2409] + v[239] * v[2411] - v[2516] * v[326] - v[2514] * v[327];
		v[2519] = i4998 + v[2499] * v[337] + v[2500] * v[339];
		v[2520] = i4998;
		v[2521] = i4999 + v[2497] * v[337] + v[2498] * v[339];
		v[2522] = i4999;
		v[2523] = v[240] * v[2413] + v[239] * v[2415] - v[2521] * v[326] - v[2519] * v[327];
		v[5048] = v[2523] * v[465] + v[2518] * v[466] + v[2513] * v[467];
		v[2524] = v[2327] + v[2331];
		v[2542] = v[1748] * v[2524];
		v[2525] = -v[2327] + v[2331];
		v[2527] = (v[189] * v[2524] + v[2323] * v[362] + v[2526] * v[5000]) / v[364];
		v[2528] = v[2325] + v[2332];
		v[2550] = v[1748] * v[2528];
		v[2529] = -v[2325] + v[2332];
		v[2546] = v[1747] * v[2529];
		v[2531] = (v[185] * v[2528] + v[2321] * v[363] + v[2530] * v[5000]) / v[364];
		v[2533] = (v[180] * v[2529] + v[2318] * v[4749] + v[2532] * v[5000]) / v[364];
		v[2534] = v[2544] + v[2546];
		v[5179] = v[2534] / v[364];
		v[2535] = v[2328] + v[2334];
		v[2540] = v[1747] * v[2535];
		v[2536] = v[2328] - v[2334];
		v[2537] = v[2548] + v[2550];
		v[5175] = v[2537] / v[364];
		v[2538] = v[1746] * v[2326] + v[2540] + v[2542];
		v[5183] = v[2538] / v[364];
		v[2541] = v[2538] - v[2542];
		v[2543] = v[2538] - v[2540];
		v[5176] = v[2543] / v[364];
		v[2545] = v[1746] * v[2525] + v[2544];
		v[2547] = v[2545] + v[2546];
		v[5177] = v[2547] / v[364];
		v[2549] = v[1746] * v[2536] + v[2548];
		v[2551] = v[2549] + v[2550];
		v[5180] = v[2551] / v[364];
		v[2552] = v[2349] + v[2353];
		v[2570] = v[1745] * v[2552];
		v[2553] = -v[2349] + v[2353];
		v[2555] = (v[208] * v[2552] + v[2345] * v[388] + v[2554] * v[5001]) / v[390];
		v[2556] = v[2347] + v[2354];
		v[2578] = v[1745] * v[2556];
		v[2557] = -v[2347] + v[2354];
		v[2574] = v[1744] * v[2557];
		v[2559] = (v[204] * v[2556] + v[2343] * v[389] + v[2558] * v[5001]) / v[390];
		v[2561] = (v[199] * v[2557] + v[2340] * v[4754] + v[2560] * v[5001]) / v[390];
		v[2562] = v[2572] + v[2574];
		v[5192] = v[2562] / v[390];
		v[2563] = v[2350] + v[2356];
		v[2568] = v[1744] * v[2563];
		v[2564] = v[2350] - v[2356];
		v[2565] = v[2576] + v[2578];
		v[5188] = v[2565] / v[390];
		v[2566] = v[1743] * v[2348] + v[2568] + v[2570];
		v[5196] = v[2566] / v[390];
		v[2569] = v[2566] - v[2570];
		v[2571] = v[2566] - v[2568];
		v[5189] = v[2571] / v[390];
		v[2573] = v[1743] * v[2553] + v[2572];
		v[2575] = v[2573] + v[2574];
		v[5190] = v[2575] / v[390];
		v[2577] = v[1743] * v[2564] + v[2576];
		v[2579] = v[2577] + v[2578];
		v[5193] = v[2579] / v[390];
		v[2580] = v[2431] + v[2435];
		v[2598] = v[1742] * v[2580];
		v[2581] = -v[2431] + v[2435];
		v[2583] = (v[2580] * v[276] + v[2427] * v[414] + v[2582] * v[5002]) / v[416];
		v[2584] = v[2429] + v[2436];
		v[2606] = v[1742] * v[2584];
		v[2585] = -v[2429] + v[2436];
		v[2602] = v[1741] * v[2585];
		v[2587] = (v[2584] * v[272] + v[2425] * v[415] + v[2586] * v[5002]) / v[416];
		v[2589] = (v[2585] * v[267] + v[2422] * v[4759] + v[2588] * v[5002]) / v[416];
		v[2590] = v[2600] + v[2602];
		v[5205] = v[2590] / v[416];
		v[2591] = v[2432] + v[2438];
		v[2596] = v[1741] * v[2591];
		v[2592] = v[2432] - v[2438];
		v[2593] = v[2604] + v[2606];
		v[5201] = v[2593] / v[416];
		v[2594] = v[1740] * v[2430] + v[2596] + v[2598];
		v[5209] = v[2594] / v[416];
		v[2597] = v[2594] - v[2598];
		v[2599] = v[2594] - v[2596];
		v[5202] = v[2599] / v[416];
		v[2601] = v[1740] * v[2581] + v[2600];
		v[2603] = v[2601] + v[2602];
		v[5203] = v[2603] / v[416];
		v[2605] = v[1740] * v[2592] + v[2604];
		v[2607] = v[2605] + v[2606];
		v[5206] = v[2607] / v[416];
		v[2608] = v[2453] + v[2457];
		v[2626] = v[1739] * v[2608];
		v[2609] = -v[2453] + v[2457];
		v[2611] = (v[2608] * v[295] + v[2449] * v[440] + v[2610] * v[5003]) / v[442];
		v[2612] = v[2451] + v[2458];
		v[2634] = v[1739] * v[2612];
		v[2613] = -v[2451] + v[2458];
		v[2630] = v[1738] * v[2613];
		v[2615] = (v[2612] * v[291] + v[2447] * v[441] + v[2614] * v[5003]) / v[442];
		v[2617] = (v[2613] * v[286] + v[2444] * v[4764] + v[2616] * v[5003]) / v[442];
		v[2618] = v[2628] + v[2630];
		v[5218] = v[2618] / v[442];
		v[2619] = v[2454] + v[2460];
		v[2624] = v[1738] * v[2619];
		v[2620] = v[2454] - v[2460];
		v[2621] = v[2632] + v[2634];
		v[5214] = v[2621] / v[442];
		v[2622] = v[1737] * v[2452] + v[2624] + v[2626];
		v[5222] = v[2622] / v[442];
		v[2625] = v[2622] - v[2626];
		v[2627] = v[2622] - v[2624];
		v[5215] = v[2627] / v[442];
		v[2629] = v[1737] * v[2609] + v[2628];
		v[2631] = v[2629] + v[2630];
		v[5216] = v[2631] / v[442];
		v[2633] = v[1737] * v[2620] + v[2632];
		v[2635] = v[2633] + v[2634];
		v[5219] = v[2635] / v[442];
		v[2636] = v[5048] / v[1475];
		v[2637] = v[2636] * v[478];
		v[2647] = v[2637] * v[467] + v[2513] * v[479];
		v[2643] = v[2637] * v[466] + v[2518] * v[479];
		v[2639] = v[2637] * v[465] + v[2523] * v[479];
		v[2640] = v[2639] * v[5004];
		v[2641] = v[2639];
		v[2644] = v[2643] * v[5005];
		v[2645] = v[2643];
		v[2648] = v[2647] * v[5006];
		v[2649] = v[2647];
		v[2650] = 0e0;
		v[2651] = 0e0;
		v[2652] = 0e0;
		v[2653] = 0e0;
		v[2654] = 0e0;
		v[2655] = 0e0;
		v[2656] = 0e0;
		v[2657] = 0e0;
		v[2658] = 0e0;
		v[2659] = 0e0;
		v[2660] = 0e0;
		v[2661] = 0e0;
		v[2662] = 0e0;
		v[2663] = 0e0;
		v[2664] = 0e0;
		v[2665] = 0e0;
		v[2666] = 0e0;
		v[2667] = 0e0;
		v[2668] = 0e0;
		v[2669] = 0e0;
		v[2670] = 0e0;
		v[2671] = 0e0;
		v[2672] = 0e0;
		v[2673] = 0e0;
		v[2674] = 0e0;
		v[2675] = 0e0;
		v[2676] = 0e0;
		v[2677] = 0e0;
		v[2678] = 0e0;
		v[2679] = 0e0;
		v[2680] = 0e0;
		v[2681] = 0e0;
		b2682 = b1307;
		if (b2682) {
			v[2683] = -(v[2645] * v[491]);
			v[2684] = v[2645] * v[489];
			v[2685] = v[2683] + v[2649] * v[490];
			v[2686] = -(v[2649] * v[489]);
			v[2687] = v[2686] + v[2641] * v[491];
			v[2688] = v[2684] - v[2641] * v[490];
			v[5010] = v[1874] * v[2685] + v[1873] * v[2687] + v[1871] * v[2688];
			v[5008] = v[1309] * v[2685] + v[1310] * v[2687] + v[1311] * v[2688];
			v[2689] = v[5008] / v[1312];
			v[5009] = v[2689] * v[5356];
			v[2699] = v[2689] * v[5007];
			v[2653] = -(v[1877] * v[2690] * v[5008]);
			v[2691] = v[2685] * v[4818] + v[1309] * v[5009];
			v[2706] = v[2691] * v[5042];
			v[2693] = v[2687] * v[4818] + v[1310] * v[5009];
			v[2710] = 2e0 * v[1322] * v[2693];
			v[2694] = v[2688] * v[4818] + v[1311] * v[5009];
			v[2707] = v[2694] * v[5039];
			v[2656] = v[2699] * v[4904] + v[1314] * v[5010];
			v[2655] = v[1318] * v[2689] * v[2697];
			v[2654] = v[1319] * v[2689] * v[4904] + v[1320] * v[5010];
			v[2680] = v[1314] * v[2699] * v[2700];
			v[2681] = v[1876] * v[2689] * v[2695] * v[2700] * v[5357];
			v[2652] = v[2688] * v[4905] + v[1871] * v[5009];
			v[2651] = v[2687] * v[4905] + v[1873] * v[5009];
			v[2650] = v[2685] * v[4905] + v[1874] * v[5009];
			v[2702] = (v[1322] * v[2691] + v[1321] * v[2693]) / 2e0;
			v[2703] = v[2706] + v[2710];
			v[2704] = v[2703] + v[2707];
			v[2705] = (v[1323] * v[2691] + v[1321] * v[2694]) / 2e0;
			v[2708] = v[2706] + v[2707];
			v[2709] = (v[1323] * v[2693] + v[1322] * v[2694]) / 2e0;
			v[2711] = v[2707] + v[2710];
			v[2659] = (v[1857] * v[2691] + v[1844] * v[2693] + 4e0 * v[2694] * v[2712]) / 2e0;
			v[2658] = (v[1864] * v[2691] + v[1844] * v[2694] + v[2693] * v[5011]) / 2e0;
			v[2657] = (v[1864] * v[2693] + v[1857] * v[2694] + v[2691] * v[5012]) / 2e0;
			v[2715] = v[2704] * v[5038];
			v[2679] = 8e0 * v[1866] * v[2704] * v[5358];
			v[2678] = v[2715] * v[5013];
			v[2716] = v[2694] + v[2702];
			v[2717] = v[2694] - v[2702];
			v[2677] = v[1772] * v[2715];
			v[2718] = -v[2693] + v[2705];
			v[2719] = v[2693] + v[2705];
			v[2676] = v[1773] * v[2715];
			v[2664] = v[1851] * v[2715] + v[1324] * v[2716];
			v[2675] = v[1774] * v[2715];
			v[2665] = (-(v[1324] * v[2708]) + v[1847] * v[2715]) / 2e0;
			v[2674] = v[2715] * v[5014];
			v[2720] = v[2691] + v[2709];
			v[2721] = -v[2691] + v[2709];
			v[2673] = v[1776] * v[2715];
			v[2667] = v[1838] * v[2715] + v[1324] * v[2718];
			v[2672] = v[1777] * v[2715];
			v[2668] = v[1834] * v[2715] + v[1324] * v[2720];
			v[2671] = v[1778] * v[2715];
			v[2669] = (-(v[1324] * v[2703]) + v[1830] * v[2715]) / 2e0;
			v[2670] = v[2715] * v[5015];
			v[2666] = v[1842] * v[2715] + v[1324] * v[2721];
			v[2663] = v[1855] * v[2715] + v[1324] * v[2719];
			v[2662] = v[1860] * v[2715] - v[1324] * v[2717];
			v[2661] = (-(v[1324] * v[2711]) + v[1865] * v[2715]) / 2e0;
			v[2660] = v[1774] * v[2716] - v[1772] * v[2717] + v[1777] * v[2718] + v[1773] * v[2719] + v[1778] * v[2720]
				+ v[1776] * v[2721] - v[2711] * v[5013] - v[2708] * v[5014] - v[2703] * v[5015];
		}
		else {
		};
		v[2722] = 0e0;
		v[2723] = 0e0;
		v[2724] = 0e0;
		v[2725] = 0e0;
		v[2726] = 0e0;
		v[2727] = 0e0;
		b2728 = b39;
		if (b2728) {
			v[5018] = -(v[1344] * v[2649]);
			v[5017] = -(v[1343] * v[2645]);
			v[5016] = -(v[1342] * v[2641]);
			v[2756] = v[1765] * v[2649];
			v[2408] = v[2408] + v[2401] * v[244];
			v[2408] = v[2408] + v[2402] * v[246];
			v[2406] = v[2406] + v[2403] * v[244];
			v[2406] = v[2406] + v[2404] * v[246];
			v[2512] = v[2512] + v[2505] * v[331];
			v[2512] = v[2512] + v[2506] * v[333];
			v[2510] = v[2510] + v[2507] * v[331];
			v[2510] = v[2510] + v[2508] * v[333];
			v[2733] = v[238] * v[2406] + v[237] * v[2408] + v[11] * v[2669] - v[2512] * v[324] - v[2510] * v[325];
			v[2669] = 0e0;
			v[2734] = v[10] * v[2668] + v[2733];
			v[2668] = 0e0;
			v[2735] = v[2734] + v[2667] * v[9];
			v[5019] = -(v[2735] * v[482]);
			v[2667] = 0e0;
			v[2412] = v[2412] + v[2397] * v[244];
			v[2412] = v[2412] + v[2398] * v[246];
			v[2410] = v[2410] + v[2399] * v[244];
			v[2410] = v[2410] + v[2400] * v[246];
			v[2517] = v[2517] + v[2501] * v[331];
			v[2517] = v[2517] + v[2502] * v[333];
			v[2515] = v[2515] + v[2503] * v[331];
			v[2515] = v[2515] + v[2504] * v[333];
			v[2740] = v[238] * v[2410] + v[237] * v[2412] + v[11] * v[2666] - v[2517] * v[324] - v[2515] * v[325];
			v[2666] = 0e0;
			v[2741] = v[10] * v[2665] + v[2740];
			v[2665] = 0e0;
			v[2742] = v[2741] + v[2664] * v[9];
			v[5020] = -(v[2742] * v[481]);
			v[2664] = 0e0;
			v[2416] = v[2416] + v[2393] * v[244];
			v[2416] = v[2416] + v[2394] * v[246];
			v[2414] = v[2414] + v[2395] * v[244];
			v[2414] = v[2414] + v[2396] * v[246];
			v[2522] = v[2522] + v[2497] * v[331];
			v[2522] = v[2522] + v[2498] * v[333];
			v[2520] = v[2520] + v[2499] * v[331];
			v[2520] = v[2520] + v[2500] * v[333];
			v[2747] = v[238] * v[2414] + v[237] * v[2416] + v[11] * v[2663] - v[2522] * v[324] - v[2520] * v[325];
			v[2663] = 0e0;
			v[2748] = v[10] * v[2662] + v[2747];
			v[2662] = 0e0;
			v[2749] = v[2748] + v[2661] * v[9];
			v[5021] = -(v[2749] * v[480]);
			v[2661] = 0e0;
			v[2750] = v[1767] * v[2641];
			v[2722] = v[2641] * v[2760];
			v[2648] = v[2648] + v[1765] * v[5016];
			v[2644] = v[2644] + v[1766] * v[5016];
			v[2641] = 0e0;
			v[2752] = v[1766] * v[2645];
			v[2753] = v[2750] + v[2752];
			v[2723] = v[2645] * v[2761];
			v[2648] = v[2648] + v[1765] * v[5017];
			v[2640] = v[2640] + v[1767] * v[5017];
			v[2645] = 0e0;
			v[2755] = v[2752] + v[2756];
			v[2757] = v[2750] + v[2756];
			v[2724] = v[2649] * v[2762];
			v[2644] = v[2644] + v[1766] * v[5018];
			v[2640] = v[2640] + v[1767] * v[5018];
			v[2649] = 0e0;
			v[2727] = v[1765] * v[2735];
			v[2726] = v[1766] * v[2742];
			v[2725] = v[1767] * v[2749];
			v[2722] = v[2722] - (2e0 * v[1767] * v[2639] + v[2755]) * v[480];
			v[2640] = v[2640] - v[1342] * v[2755] + v[2749] * v[2760] + v[1767] * (v[5019] + v[5020]);
			v[2723] = v[2723] - (2e0 * v[1766] * v[2643] + v[2757]) * v[481];
			v[2644] = v[2644] - v[1343] * v[2757] + v[2742] * v[2761] + v[1766] * (v[5019] + v[5021]);
			v[2724] = v[2724] - (2e0 * v[1765] * v[2647] + v[2753]) * v[482];
			v[2648] = v[2648] - v[1344] * v[2753] + v[2735] * v[2762] + v[1765] * (v[5020] + v[5021]);
		}
		else {
		};
		v[2767] = v[2329] * v[2763] + v[2333] * v[2764] + v[2524] * v[2765] + v[2528] * v[2766] + v[2527] * v[373]
			+ v[2531] * v[375] + v[2358] * v[3976] + v[2362] * v[4750] + v[2367] * v[4752] + v[8005 + i1663];
		v[2772] = v[2210] * v[2357] + v[2209] * v[2369] + v[2329] * v[2768] + v[2330] * v[2769] + v[2529] * v[2770]
			+ v[2535] * v[2771] + v[2527] * v[367] + v[2533] * v[375] + v[2363] * v[3978] + v[8053 + i1663];
		v[2777] = v[2212] * v[2361] + v[2211] * v[2366] + v[2329] * v[2773] + v[2326] * v[2774] + v[2525] * v[2775]
			+ v[2536] * v[2776] + v[2531] * v[367] + v[2533] * v[373] + v[2371] * v[3980] + v[8101 + i1663];
		v[2782] = v[2351] * v[2778] + v[2355] * v[2779] + v[2552] * v[2780] + v[2556] * v[2781] + v[2373] * v[3982]
			+ v[2555] * v[399] + v[2559] * v[401] + v[2377] * v[4755] + v[2382] * v[4757] + v[8149 + i1663];
		v[2787] = v[2216] * v[2372] + v[2215] * v[2384] + v[2351] * v[2783] + v[2352] * v[2784] + v[2557] * v[2785]
			+ v[2563] * v[2786] + v[2555] * v[393] + v[2378] * v[3984] + v[2561] * v[401] + v[8197 + i1663];
		v[2792] = v[2218] * v[2376] + v[2217] * v[2381] + v[2351] * v[2788] + v[2348] * v[2789] + v[2553] * v[2790]
			+ v[2564] * v[2791] + v[2559] * v[393] + v[2386] * v[3986] + v[2561] * v[399] + v[8245 + i1663];
		v[2797] = v[2433] * v[2793] + v[2437] * v[2794] + v[2580] * v[2795] + v[2584] * v[2796] + v[2462] * v[3988]
			+ v[2583] * v[425] + v[2587] * v[427] + v[2466] * v[4760] + v[2471] * v[4762] + v[8293 + i1663];
		v[2802] = v[2222] * v[2461] + v[2221] * v[2473] + v[2433] * v[2798] + v[2434] * v[2799] + v[2585] * v[2800]
			+ v[2591] * v[2801] + v[2467] * v[3990] + v[2583] * v[419] + v[2589] * v[427] + v[8341 + i1663];
		v[2807] = v[2224] * v[2465] + v[2223] * v[2470] + v[2433] * v[2803] + v[2430] * v[2804] + v[2581] * v[2805]
			+ v[2592] * v[2806] + v[2475] * v[3992] + v[2587] * v[419] + v[2589] * v[425] + v[8389 + i1663];
		v[2812] = v[2455] * v[2808] + v[2459] * v[2809] + v[2608] * v[2810] + v[2612] * v[2811] + v[2477] * v[3994]
			+ v[2611] * v[451] + v[2615] * v[453] + v[2481] * v[4765] + v[2486] * v[4767] + v[8437 + i1663];
		v[2817] = v[2228] * v[2476] + v[2227] * v[2488] + v[2455] * v[2813] + v[2456] * v[2814] + v[2613] * v[2815]
			+ v[2619] * v[2816] + v[2482] * v[3996] + v[2611] * v[445] + v[2617] * v[453] + v[8485 + i1663];
		v[2822] = v[2230] * v[2480] + v[2229] * v[2485] + v[2455] * v[2818] + v[2452] * v[2819] + v[2609] * v[2820]
			+ v[2620] * v[2821] + v[2490] * v[3998] + v[2615] * v[445] + v[2617] * v[451] + v[8533 + i1663];
		v[2823] = v[18] * (v[239] * ((v[1665] * v[2393] + v[1666] * v[2397] + v[1667] * v[2401]) * v[508] + (v[1665] * v[2394]
			+ v[1666] * v[2398] + v[1667] * v[2402]) * v[509]) + v[240] * ((v[1665] * v[2395] + v[1666] * v[2399] + v[1667] * v[2403]
				) * v[508] + (v[1665] * v[2396] + v[1666] * v[2400] + v[1667] * v[2404]) * v[509]));
		v[2824] = (v[2405] - v[2407]) * v[2826] + (v[2409] - v[2411]) * v[2827] + (v[2413] - v[2415]) * v[2828];
		v[2825] = v[18] * (v[326] * ((v[1665] * v[2497] + v[1666] * v[2501] + v[1667] * v[2505]) * v[536] + (v[1665] * v[2498]
			+ v[1666] * v[2502] + v[1667] * v[2506]) * v[537]) + v[327] * ((v[1665] * v[2499] + v[1666] * v[2503] + v[1667] * v[2507]
				) * v[536] + (v[1665] * v[2500] + v[1666] * v[2504] + v[1667] * v[2508]) * v[537]));
		v[2829] = (v[2509] - v[2511]) * v[2826] + (v[2514] - v[2516]) * v[2827] + (v[2519] - v[2521]) * v[2828];
		v[2830] = -(i4985 * v[1566]) - i4983 * v[1570] - i4981 * v[1574] - i4984 * v[1590] - i4982 * v[1594] - i4980 * v[1598]
			- i4999 * v[1614] - i4997 * v[1618] - i4995 * v[1622] - i4998 * v[1638] - i4996 * v[1642] - i4994 * v[1646]
			- v[1658] * v[2277] - v[1654] * v[2278] - v[1650] * v[2279] - v[1634] * v[2280] - v[1630] * v[2281] - v[1626] * v[2282]
			- v[1610] * v[2283] - v[1606] * v[2284] - v[1602] * v[2285] - v[1586] * v[2286] - v[1582] * v[2287] - v[1578] * v[2288];
		v[2845] = v[2830];
		v[2831] = -(i4985 * v[1567]) - i4983 * v[1571] - i4981 * v[1575] - i4984 * v[1591] - i4982 * v[1595] - i4980 * v[1599]
			- i4999 * v[1615] - i4997 * v[1619] - i4995 * v[1623] - i4998 * v[1639] - i4996 * v[1643] - i4994 * v[1647]
			- v[1659] * v[2277] - v[1655] * v[2278] - v[1651] * v[2279] - v[1635] * v[2280] - v[1631] * v[2281] - v[1627] * v[2282]
			- v[1611] * v[2283] - v[1607] * v[2284] - v[1603] * v[2285] - v[1587] * v[2286] - v[1583] * v[2287] - v[1579] * v[2288];
		v[2843] = v[2831];
		v[2832] = -(i4985 * v[1568]) - i4983 * v[1572] - i4981 * v[1576] - i4984 * v[1592] - i4982 * v[1596] - i4980 * v[1600]
			- i4999 * v[1616] - i4997 * v[1620] - i4995 * v[1624] - i4998 * v[1640] - i4996 * v[1644] - i4994 * v[1648]
			- v[1660] * v[2277] - v[1656] * v[2278] - v[1652] * v[2279] - v[1636] * v[2280] - v[1632] * v[2281] - v[1628] * v[2282]
			- v[1612] * v[2283] - v[1608] * v[2284] - v[1604] * v[2285] - v[1588] * v[2286] - v[1584] * v[2287] - v[1580] * v[2288];
		v[2841] = v[2832];
		v[2833] = 0e0;
		v[2834] = 0e0;
		v[2835] = 0e0;
		v[2836] = 0e0;
		v[2837] = 0e0;
		v[2838] = 0e0;
		b2839 = b38;
		if (b2839) {
			b2840 = b1496;
			if (b2840) {
				v[2838] = v[2832];
				v[2832] = 0e0;
				v[2837] = v[2831];
				v[2831] = 0e0;
				v[2836] = v[2830];
				v[2830] = 0e0;
			}
			else {
				v[2853] = v[2845] * v[4894];
				v[2851] = v[2843] * v[4894];
				v[2849] = v[2841] * v[4894];
				v[2832] = 0e0;
				v[2831] = 0e0;
				v[5023] = (v[16] * (v[1510] * v[2841] + v[1509] * v[2843] + v[1508] * v[2845])) / sqrt(v[2854]);
				v[2830] = 0e0;
				v[5022] = (v[1506] * (v[1494] * v[2849] + v[1493] * v[2851] + v[1492] * v[2853])) / sqrt(v[2850]);
				v[2838] = v[1507] * v[2849] + v[1494] * v[5022];
				v[2837] = v[1507] * v[2851] + v[1493] * v[5022];
				v[2836] = v[1507] * v[2853] + v[1492] * v[5022];
				v[2835] = v[1485] * v[5023];
				v[2834] = v[1484] * v[5023];
				v[2833] = v[1483] * v[5023];
			};
		}
		else {
			b2856 = b1530;
			if (b2856) {
				v[2838] = v[2841];
				v[2832] = 0e0;
				v[2837] = v[2843];
				v[2831] = 0e0;
				v[2836] = v[2845];
				v[2830] = 0e0;
			}
			else {
				v[2862] = v[1541] * v[16] * v[2845];
				v[2860] = v[1541] * v[16] * v[2843];
				v[2859] = v[1541] * v[16] * v[2841];
				v[5025] = (v[16] * (v[1540] * v[2841] + v[1539] * v[2843] + v[1538] * v[2845])) / sqrt(v[2854]);
				v[5024] = (v[1536] * (v[1494] * v[2859] + v[1493] * v[2860] + v[1492] * v[2862])) / sqrt(v[2850]);
				v[2838] = v[1537] * v[2859] + v[1494] * v[5024];
				v[2837] = v[1537] * v[2860] + v[1493] * v[5024];
				v[2836] = v[1537] * v[2862] + v[1492] * v[5024];
				v[2835] = v[1485] * v[5025];
				v[2834] = v[1484] * v[5025];
				v[2833] = v[1483] * v[5025];
			};
		};
		v[5026] = v[17] * v[2833];
		v[3212] = v[1383] * v[5026];
		v[3206] = v[480] * v[5026];
		v[5027] = v[17] * v[2834];
		v[3112] = v[1413] * v[5027];
		v[3109] = v[481] * v[5027];
		v[5028] = v[17] * v[2835];
		v[3021] = v[2835] * v[4773];
		v[5092] = -v[3212] - v[3021] * v[480];
		v[5091] = -v[3112] - v[3021] * v[481];
		v[5031] = v[3212] + (v[3021] + v[3109]) * v[480];
		v[5030] = v[3112] + (v[3021] + v[3206]) * v[481];
		v[5029] = (v[3109] + v[3206]) * v[482] + v[1455] * v[5028];
		v[2868] = -(v[1500] * v[2277]) - v[18] * (v[1667] * v[2822] + v[2838] * v[464]);
		v[2869] = -(v[1500] * v[2278]) - v[18] * (v[1667] * v[2817] + v[2838] * v[454]);
		v[2870] = -(v[1500] * v[2279]) - v[18] * (v[1667] * v[2812] + v[2838] * v[444]);
		v[2871] = -(i4994 * v[1500]) - v[18] * (v[1667] * v[2496] + v[2838] * v[360]);
		v[2872] = -(i4996 * v[1500]) - v[18] * (v[1667] * v[2495] + v[2838] * v[359]);
		v[2873] = -(i4998 * v[1500]) - v[18] * (v[1667] * v[2494] + v[2838] * v[358]);
		v[2874] = -(v[1500] * v[2280]) - v[18] * (v[1667] * v[2807] + v[2838] * v[438]);
		v[2875] = -(v[1500] * v[2281]) - v[18] * (v[1667] * v[2802] + v[2838] * v[428]);
		v[2876] = -(v[1500] * v[2282]) - v[18] * (v[1667] * v[2797] + v[2838] * v[418]);
		v[2877] = -(i4995 * v[1500]) - v[18] * (v[1667] * v[2493] + v[2838] * v[357]);
		v[2878] = -(i4997 * v[1500]) - v[18] * (v[1667] * v[2492] + v[2838] * v[356]);
		v[2879] = -(i4999 * v[1500]) - v[18] * (v[1667] * v[2491] + v[2838] * v[355]);
		v[2880] = -(v[1500] * v[2283]) - v[18] * (v[1667] * v[2792] + v[2838] * v[412]);
		v[2881] = -(v[1500] * v[2284]) - v[18] * (v[1667] * v[2787] + v[2838] * v[402]);
		v[2882] = -(v[1500] * v[2285]) - v[18] * (v[1667] * v[2782] + v[2838] * v[392]);
		v[2883] = -(i4980 * v[1500]) - v[18] * (v[1667] * v[2392] + v[2838] * v[354]);
		v[2884] = -(i4982 * v[1500]) - v[18] * (v[1667] * v[2391] + v[2838] * v[353]);
		v[2885] = -(i4984 * v[1500]) - v[18] * (v[1667] * v[2390] + v[2838] * v[352]);
		v[2886] = -(v[1500] * v[2286]) - v[18] * (v[1667] * v[2777] + v[2838] * v[386]);
		v[2887] = -(v[1500] * v[2287]) - v[18] * (v[1667] * v[2772] + v[2838] * v[376]);
		v[2888] = -(v[1500] * v[2288]) - v[18] * (v[1667] * v[2767] + v[2838] * v[366]);
		v[2889] = -(i4981 * v[1500]) - v[18] * (v[1667] * v[2389] + v[2838] * v[351]);
		v[2890] = -(i4983 * v[1500]) - v[18] * (v[1667] * v[2388] + v[2838] * v[350]);
		v[2891] = -(i4985 * v[1500]) - v[18] * (v[1667] * v[2387] + v[2838] * v[349]);
		v[2917] = -(v[1499] * v[2277]) - v[18] * (v[1666] * v[2822] + v[2837] * v[464]);
		v[2918] = -(v[1499] * v[2278]) - v[18] * (v[1666] * v[2817] + v[2837] * v[454]);
		v[2919] = -(v[1499] * v[2279]) - v[18] * (v[1666] * v[2812] + v[2837] * v[444]);
		v[2920] = -(i4994 * v[1499]) - v[18] * (v[1666] * v[2496] + v[2837] * v[360]);
		v[2921] = -(i4996 * v[1499]) - v[18] * (v[1666] * v[2495] + v[2837] * v[359]);
		v[2922] = -(i4998 * v[1499]) - v[18] * (v[1666] * v[2494] + v[2837] * v[358]);
		v[2923] = -(v[1499] * v[2280]) - v[18] * (v[1666] * v[2807] + v[2837] * v[438]);
		v[2924] = -(v[1499] * v[2281]) - v[18] * (v[1666] * v[2802] + v[2837] * v[428]);
		v[2925] = -(v[1499] * v[2282]) - v[18] * (v[1666] * v[2797] + v[2837] * v[418]);
		v[2926] = -(i4995 * v[1499]) - v[18] * (v[1666] * v[2493] + v[2837] * v[357]);
		v[2927] = -(i4997 * v[1499]) - v[18] * (v[1666] * v[2492] + v[2837] * v[356]);
		v[2928] = -(i4999 * v[1499]) - v[18] * (v[1666] * v[2491] + v[2837] * v[355]);
		v[2929] = -(v[1499] * v[2283]) - v[18] * (v[1666] * v[2792] + v[2837] * v[412]);
		v[2930] = -(v[1499] * v[2284]) - v[18] * (v[1666] * v[2787] + v[2837] * v[402]);
		v[2931] = -(v[1499] * v[2285]) - v[18] * (v[1666] * v[2782] + v[2837] * v[392]);
		v[2932] = -(i4980 * v[1499]) - v[18] * (v[1666] * v[2392] + v[2837] * v[354]);
		v[2933] = -(i4982 * v[1499]) - v[18] * (v[1666] * v[2391] + v[2837] * v[353]);
		v[2934] = -(i4984 * v[1499]) - v[18] * (v[1666] * v[2390] + v[2837] * v[352]);
		v[2935] = -(v[1499] * v[2286]) - v[18] * (v[1666] * v[2777] + v[2837] * v[386]);
		v[2936] = -(v[1499] * v[2287]) - v[18] * (v[1666] * v[2772] + v[2837] * v[376]);
		v[2937] = -(v[1499] * v[2288]) - v[18] * (v[1666] * v[2767] + v[2837] * v[366]);
		v[2938] = -(i4981 * v[1499]) - v[18] * (v[1666] * v[2389] + v[2837] * v[351]);
		v[2939] = -(i4983 * v[1499]) - v[18] * (v[1666] * v[2388] + v[2837] * v[350]);
		v[2940] = -(i4985 * v[1499]) - v[18] * (v[1666] * v[2387] + v[2837] * v[349]);
		v[2966] = -(v[1498] * v[2277]) - v[18] * (v[1665] * v[2822] + v[2836] * v[464]);
		v[2967] = -(v[1498] * v[2278]) - v[18] * (v[1665] * v[2817] + v[2836] * v[454]);
		v[2968] = -(v[1498] * v[2279]) - v[18] * (v[1665] * v[2812] + v[2836] * v[444]);
		v[2969] = -(i4994 * v[1498]) - v[18] * (v[1665] * v[2496] + v[2836] * v[360]);
		v[2970] = -(i4996 * v[1498]) - v[18] * (v[1665] * v[2495] + v[2836] * v[359]);
		v[2971] = -(i4998 * v[1498]) - v[18] * (v[1665] * v[2494] + v[2836] * v[358]);
		v[2972] = -(v[1498] * v[2280]) - v[18] * (v[1665] * v[2807] + v[2836] * v[438]);
		v[2973] = -(v[1498] * v[2281]) - v[18] * (v[1665] * v[2802] + v[2836] * v[428]);
		v[2974] = -(v[1498] * v[2282]) - v[18] * (v[1665] * v[2797] + v[2836] * v[418]);
		v[2975] = -(i4995 * v[1498]) - v[18] * (v[1665] * v[2493] + v[2836] * v[357]);
		v[2976] = -(i4997 * v[1498]) - v[18] * (v[1665] * v[2492] + v[2836] * v[356]);
		v[2977] = -(i4999 * v[1498]) - v[18] * (v[1665] * v[2491] + v[2836] * v[355]);
		v[2978] = -(v[1498] * v[2283]) - v[18] * (v[1665] * v[2792] + v[2836] * v[412]);
		v[2979] = -(v[1498] * v[2284]) - v[18] * (v[1665] * v[2787] + v[2836] * v[402]);
		v[2980] = -(v[1498] * v[2285]) - v[18] * (v[1665] * v[2782] + v[2836] * v[392]);
		v[2981] = -(i4980 * v[1498]) - v[18] * (v[1665] * v[2392] + v[2836] * v[354]);
		v[2982] = -(i4982 * v[1498]) - v[18] * (v[1665] * v[2391] + v[2836] * v[353]);
		v[2983] = -(i4984 * v[1498]) - v[18] * (v[1665] * v[2390] + v[2836] * v[352]);
		v[2984] = -(v[1498] * v[2286]) - v[18] * (v[1665] * v[2777] + v[2836] * v[386]);
		v[2985] = -(v[1498] * v[2287]) - v[18] * (v[1665] * v[2772] + v[2836] * v[376]);
		v[2986] = -(v[1498] * v[2288]) - v[18] * (v[1665] * v[2767] + v[2836] * v[366]);
		v[2987] = -(i4981 * v[1498]) - v[18] * (v[1665] * v[2389] + v[2836] * v[351]);
		v[2988] = -(i4983 * v[1498]) - v[18] * (v[1665] * v[2388] + v[2836] * v[350]);
		v[2989] = -(i4985 * v[1498]) - v[18] * (v[1665] * v[2387] + v[2836] * v[349]);
		v[3014] = v[14] * v[2838];
		v[3015] = v[14] * v[2837];
		v[3016] = v[14] * v[2836];
		v[3055] = v[4388] * v[5028];
		v[2648] = v[2648] + v[2835] * v[4389];
		v[2644] = v[2644] + v[4395] * v[5028];
		v[2640] = v[2640] + v[4401] * v[5028];
		v[2648] = v[2648] + v[4387] * v[5027];
		v[3152] = v[4847] * v[5027];
		v[2644] = v[2644] + v[2834] * v[4394];
		v[2640] = v[2640] + v[4400] * v[5027];
		v[3203] = v[376] * v[5029];
		v[3205] = v[366] * v[5029];
		v[3207] = v[376] * v[5030];
		v[3208] = v[366] * v[5031];
		v[3209] = v[386] * v[5029];
		v[3210] = v[386] * v[5030];
		v[3211] = v[366] * v[5030];
		v[3213] = v[386] * v[5031];
		v[3214] = v[376] * v[5031];
		v[3215] = v[402] * v[5029];
		v[3216] = v[392] * v[5029];
		v[3217] = v[402] * v[5030];
		v[3218] = v[392] * v[5031];
		v[3219] = v[412] * v[5029];
		v[3220] = v[412] * v[5030];
		v[3221] = v[392] * v[5030];
		v[3222] = v[412] * v[5031];
		v[3223] = v[402] * v[5031];
		v[3224] = v[428] * v[5029];
		v[3225] = v[418] * v[5029];
		v[3226] = v[428] * v[5030];
		v[3227] = v[418] * v[5031];
		v[3228] = v[438] * v[5029];
		v[3229] = v[438] * v[5030];
		v[3230] = v[418] * v[5030];
		v[3231] = v[438] * v[5031];
		v[3232] = v[428] * v[5031];
		v[3233] = v[454] * v[5029];
		v[3234] = v[444] * v[5029];
		v[3235] = v[454] * v[5030];
		v[3236] = v[444] * v[5031];
		v[3237] = v[464] * v[5029];
		v[3238] = v[464] * v[5030];
		v[3239] = v[444] * v[5030];
		v[3240] = v[464] * v[5031];
		v[3241] = v[454] * v[5031];
		v[3242] = v[3109] + v[3206];
		v[5089] = v[3242] * v[482];
		v[2648] = v[2648] + v[4386] * v[5026];
		v[3246] = v[17] * (v[2834] * v[349] + v[2833] * v[350]);
		v[3247] = v[17] * (v[2834] * v[352] + v[2833] * v[353]);
		v[3248] = v[17] * (v[2834] * v[355] + v[2833] * v[356]);
		v[3249] = v[17] * (v[2834] * v[358] + v[2833] * v[359]);
		v[5047] = v[239] * v[3246] + v[240] * v[3247] - v[3248] * v[326] - v[3249] * v[327];
		v[2644] = v[2644] + v[4391] * v[5026];
		v[3250] = v[4845] * v[5026];
		v[2640] = v[2640] + v[2833] * v[4399];
		v[3313] = v[3203] * v[985] + v[3205] * v[986];
		v[5162] = v[3313] + v[3209] * v[981];
		v[3314] = v[3207] + v[3208];
		v[3315] = v[3207] + v[3209];
		v[3316] = v[3208] + v[3209];
		v[3317] = v[5162] * v[962];
		v[3318] = v[3210] * v[981] + v[3211] * v[986];
		v[5163] = v[3318] + v[3207] * v[985];
		v[3319] = v[5163] * v[962];
		v[3320] = v[3213] * v[981] + v[3214] * v[985];
		v[5164] = v[3320] + v[3208] * v[986];
		v[3321] = v[3205] - v[3213] - v[3318] / 2e0 + v[3214] * v[982] + v[3203] * v[983] + v[3316] * v[984];
		v[3322] = -v[3211] + v[3214] - v[3313] / 2e0 + v[3213] * v[982] + v[3314] * v[983] + v[3210] * v[984];
		v[3323] = v[5164] * v[962];
		v[3324] = -v[3203] + v[3210] - v[3320] / 2e0 + v[3315] * v[982] + v[3205] * v[983] + v[3211] * v[984];
		v[3325] = v[3215] * v[979] + v[3216] * v[980];
		v[5149] = v[3325] + v[3219] * v[975];
		v[3326] = v[3217] + v[3218];
		v[3327] = v[3217] + v[3219];
		v[3328] = v[3218] + v[3219];
		v[3329] = v[5149] * v[931];
		v[3330] = v[3220] * v[975] + v[3221] * v[980];
		v[5150] = v[3330] + v[3217] * v[979];
		v[3331] = v[5150] * v[931];
		v[3332] = v[3222] * v[975] + v[3223] * v[979];
		v[5151] = v[3332] + v[3218] * v[980];
		v[3333] = v[3216] - v[3222] - v[3330] / 2e0 + v[3223] * v[976] + v[3215] * v[977] + v[3328] * v[978];
		v[3334] = -v[3221] + v[3223] - v[3325] / 2e0 + v[3222] * v[976] + v[3326] * v[977] + v[3220] * v[978];
		v[3335] = v[5151] * v[931];
		v[3336] = -v[3215] + v[3220] - v[3332] / 2e0 + v[3327] * v[976] + v[3216] * v[977] + v[3221] * v[978];
		v[3337] = v[3224] * v[973] + v[3225] * v[974];
		v[5136] = v[3337] + v[3228] * v[969];
		v[3338] = v[3226] + v[3227];
		v[3339] = v[3226] + v[3228];
		v[3340] = v[3227] + v[3228];
		v[3341] = v[5136] * v[780];
		v[3342] = v[3229] * v[969] + v[3230] * v[974];
		v[5137] = v[3342] + v[3226] * v[973];
		v[3343] = v[5137] * v[780];
		v[3344] = v[3231] * v[969] + v[3232] * v[973];
		v[5138] = v[3344] + v[3227] * v[974];
		v[3345] = v[3225] - v[3231] - v[3342] / 2e0 + v[3232] * v[970] + v[3224] * v[971] + v[3340] * v[972];
		v[3346] = -v[3230] + v[3232] - v[3337] / 2e0 + v[3231] * v[970] + v[3338] * v[971] + v[3229] * v[972];
		v[3347] = v[5138] * v[780];
		v[3348] = -v[3224] + v[3229] - v[3344] / 2e0 + v[3339] * v[970] + v[3225] * v[971] + v[3230] * v[972];
		v[3349] = v[3233] * v[967] + v[3234] * v[968];
		v[5123] = v[3349] + v[3237] * v[963];
		v[3350] = v[3235] + v[3236];
		v[3351] = v[3235] + v[3237];
		v[3352] = v[3236] + v[3237];
		v[3353] = v[5123] * v[749];
		v[3354] = v[3238] * v[963] + v[3239] * v[968];
		v[5124] = v[3354] + v[3235] * v[967];
		v[3355] = v[5124] * v[749];
		v[3356] = v[3240] * v[963] + v[3241] * v[967];
		v[5125] = v[3356] + v[3236] * v[968];
		v[3357] = v[3234] - v[3240] - v[3354] / 2e0 + v[3241] * v[964] + v[3233] * v[965] + v[3352] * v[966];
		v[3358] = -v[3239] + v[3241] - v[3349] / 2e0 + v[3240] * v[964] + v[3350] * v[965] + v[3238] * v[966];
		v[3359] = v[5125] * v[749];
		v[3360] = -v[3233] + v[3238] - v[3356] / 2e0 + v[3351] * v[964] + v[3234] * v[965] + v[3239] * v[966];
		v[3361] = 0e0;
		v[3362] = 0e0;
		v[3363] = 0e0;
		v[3364] = 0e0;
		v[3365] = 0e0;
		v[3366] = 0e0;
		v[3367] = 0e0;
		v[3368] = 0e0;
		v[3369] = 0e0;
		b3370 = b39;
		if (b3370) {
			v[5034] = v[3016] * v[480];
			v[5033] = v[3015] * v[481];
			v[5035] = v[5033] + v[5034];
			v[5032] = v[3014] * v[482];
			v[5037] = v[5032] + v[5033];
			v[5036] = v[5032] + v[5034];
			v[2727] = v[2727] + v[1344] * v[3014];
			v[2726] = v[2726] + v[1343] * v[3015];
			v[2722] = v[2722] + v[1783] * v[3016] - v[480] * v[5037];
			v[2723] = v[2723] + v[1785] * v[3015] - v[481] * v[5036];
			v[2724] = v[2724] + v[1787] * v[3014] - v[482] * v[5035];
			v[2725] = v[2725] + v[1342] * v[3016];
			v[2648] = v[2648] + v[3014] * v[4900] - v[1344] * v[5035];
			v[2644] = v[2644] + v[3015] * v[4901] - v[1343] * v[5036];
			v[2640] = v[2640] - v[3016] * v[4902] - v[1342] * v[5037];
			v[3361] = v[2722] * v[9];
			v[3362] = v[10] * v[2722];
			v[3363] = v[11] * v[2722];
			v[3371] = -(v[2722] * v[325]);
			v[3372] = -(v[2722] * v[324]);
			v[3373] = v[333] * v[3371];
			v[3374] = v[331] * v[3371];
			v[3375] = v[333] * v[3372];
			v[3376] = v[331] * v[3372];
			v[3377] = v[238] * v[2722];
			v[3378] = v[237] * v[2722];
			v[3379] = v[246] * v[3377];
			v[3380] = v[244] * v[3377];
			v[3381] = v[246] * v[3378];
			v[3382] = v[244] * v[3378];
			v[3364] = v[2723] * v[9];
			v[3365] = v[10] * v[2723];
			v[3366] = v[11] * v[2723];
			v[3383] = -(v[2723] * v[325]);
			v[3384] = -(v[2723] * v[324]);
			v[3385] = v[333] * v[3383];
			v[3386] = v[331] * v[3383];
			v[3387] = v[333] * v[3384];
			v[3388] = v[331] * v[3384];
			v[3389] = v[238] * v[2723];
			v[3390] = v[237] * v[2723];
			v[3391] = v[246] * v[3389];
			v[3392] = v[244] * v[3389];
			v[3393] = v[246] * v[3390];
			v[3394] = v[244] * v[3390];
			v[3367] = v[2724] * v[9];
			v[3368] = v[10] * v[2724];
			v[3369] = v[11] * v[2724];
			v[3395] = -(v[2724] * v[325]);
			v[3396] = -(v[2724] * v[324]);
			v[3397] = v[333] * v[3395];
			v[3398] = v[331] * v[3395];
			v[3399] = v[333] * v[3396];
			v[3400] = v[331] * v[3396];
			v[3401] = v[238] * v[2724];
			v[3402] = v[237] * v[2724];
			v[3403] = v[246] * v[3401];
			v[3404] = v[244] * v[3401];
			v[3405] = v[246] * v[3402];
			v[3406] = v[244] * v[3402];
			v[3250] = -v[2725] + v[3250];
			v[3152] = -v[2726] + v[3152];
			v[3055] = -v[2727] + v[3055];
		}
		else {
			v[3382] = 0e0;
			v[3381] = 0e0;
			v[3394] = 0e0;
			v[3393] = 0e0;
			v[3406] = 0e0;
			v[3405] = 0e0;
			v[3380] = 0e0;
			v[3379] = 0e0;
			v[3392] = 0e0;
			v[3391] = 0e0;
			v[3404] = 0e0;
			v[3403] = 0e0;
			v[3378] = 0e0;
			v[3390] = 0e0;
			v[3402] = 0e0;
			v[3377] = 0e0;
			v[3389] = 0e0;
			v[3401] = 0e0;
			v[3376] = 0e0;
			v[3375] = 0e0;
			v[3388] = 0e0;
			v[3387] = 0e0;
			v[3400] = 0e0;
			v[3399] = 0e0;
			v[3374] = 0e0;
			v[3373] = 0e0;
			v[3386] = 0e0;
			v[3385] = 0e0;
			v[3398] = 0e0;
			v[3397] = 0e0;
			v[3372] = 0e0;
			v[3384] = 0e0;
			v[3396] = 0e0;
			v[3371] = 0e0;
			v[3383] = 0e0;
			v[3395] = 0e0;
		};
		b3407 = b1307;
		if (b3407) {
			v[2660] = v[2660] + (v[1830] * v[3369]) / 2e0;
			v[2670] = v[2670] + v[3369] * v[4820];
			v[2660] = v[2660] + v[1834] * v[3368];
			v[2671] = v[2671] + v[1324] * v[3368];
			v[2660] = v[2660] + v[1838] * v[3367];
			v[2672] = v[2672] + v[1324] * v[3367];
			v[2660] = v[2660] + v[1842] * v[3366];
			v[2673] = v[2673] + v[1324] * v[3366];
			v[2660] = v[2660] + (v[1847] * v[3365]) / 2e0;
			v[2674] = v[2674] + v[3365] * v[4820];
			v[2660] = v[2660] + v[1851] * v[3364];
			v[2675] = v[2675] + v[1324] * v[3364];
			v[2660] = v[2660] + v[1855] * v[3363];
			v[2676] = v[2676] + v[1324] * v[3363];
			v[2660] = v[2660] + v[1860] * v[3362];
			v[2677] = v[2677] + v[1324] * v[3362];
			v[2660] = v[2660] + (v[1865] * v[3361]) / 2e0;
			v[2678] = v[2678] + v[3361] * v[4820];
			v[2679] = v[2679] + v[2660] * v[5038];
			v[5043] = -v[2674] + v[2679];
			v[2658] = v[2658] - v[2672];
			v[3417] = v[2672] + v[2676];
			v[2658] = v[2658] + v[2676];
			v[2657] = v[2657] + v[2671] + (v[1323] * v[3417]) / 2e0;
			v[3419] = v[2671] + v[2673];
			v[2657] = v[2657] - v[2673];
			v[2659] = v[2659] + v[2675] + v[3419] * v[4819] + v[3417] * v[4903] + v[5039] * (-v[2678] + v[5043]);
			v[2659] = v[2659] - v[2677];
			v[5040] = v[1311] * v[2659];
			v[3421] = v[2675] + v[2677];
			v[2656] = v[2656] + v[1314] * v[5040];
			v[2654] = v[2654] + v[1320] * v[5040];
			v[2652] = v[2652] + v[2659] * v[4818];
			v[2658] = v[2658] + (-4e0 * v[1322] * (v[2670] + v[2678] - v[2679]) + v[1323] * v[3419] + v[1321] * v[3421]) / 2e0;
			v[5041] = v[1310] * v[2658];
			v[2656] = v[2656] + v[1314] * v[5041];
			v[2654] = v[2654] + v[1320] * v[5041];
			v[2651] = v[2651] + v[2658] * v[4818];
			v[2657] = v[2657] + v[3421] * v[4819] + v[5042] * (-v[2670] + v[5043]);
			v[5044] = v[1309] * v[2657];
			v[2656] = v[2656] + v[1314] * v[5044];
			v[2654] = v[2654] + v[1320] * v[5044];
			v[2650] = v[2650] + v[2657] * v[4818];
			v[2655] = v[2655] + v[1319] * v[2656];
			v[2653] = v[2653] + v[2655];
			v[2680] = v[2680] + 2e0 * v[2654] * v[2695];
			v[2681] = v[2681] + (v[2680] * v[2696]) / 2e0;
			v[2653] = v[2653] + v[2681];
			v[5045] = v[2653] / v[1312];
			v[2652] = v[2652] + v[1311] * v[5045];
			v[2651] = v[2651] + v[1310] * v[5045];
			v[2650] = v[2650] + v[1309] * v[5045];
			v[2644] = v[2644] + v[2652] * v[489];
			v[2640] = v[2640] - v[2652] * v[490];
			v[2648] = v[2648] - v[2651] * v[489];
			v[2640] = v[2640] + v[2651] * v[491];
			v[2648] = v[2648] + v[2650] * v[490];
			v[2644] = v[2644] - v[2650] * v[491];
		}
		else {
		};
		v[3610] = (v[1282] * v[2966] + v[1281] * v[2967] + v[1280] * v[2968] + v[1279] * v[2969] + v[1278] * v[2970]
			+ v[1277] * v[2971] + v[1276] * v[2972] + v[1275] * v[2973] + v[1274] * v[2974] + v[1273] * v[2975] + v[1272] * v[2976]
			+ v[1271] * v[2977] + v[1270] * v[2978] + v[1269] * v[2979] + v[1268] * v[2980] + v[1267] * v[2981] + v[1266] * v[2982]
			+ v[1265] * v[2983] + v[1264] * v[2984] + v[1263] * v[2985] + v[1262] * v[2986] + v[1261] * v[2987] + v[1260] * v[2988]
			+ v[1259] * v[2989]) / 2e0;
		v[3428] = -(v[1306] * v[2966]) - v[1305] * v[2967] - v[1304] * v[2968] - v[1303] * v[2969] - v[1302] * v[2970]
			- v[1301] * v[2971] - v[1300] * v[2972] - v[1299] * v[2973] - v[1298] * v[2974] - v[1297] * v[2975] - v[1296] * v[2976]
			- v[1295] * v[2977] - v[1294] * v[2978] - v[1293] * v[2979] - v[1292] * v[2980] - v[1291] * v[2981] - v[1290] * v[2982]
			- v[1289] * v[2983] - v[1288] * v[2984] - v[1287] * v[2985] - v[1286] * v[2986] - v[1285] * v[2987] - v[1284] * v[2988]
			- v[1283] * v[2989];
		v[5086] = v[326] * v[3428];
		v[5085] = v[327] * v[3428];
		v[3744] = (-(v[1234] * v[2966]) - v[1233] * v[2967] - v[1232] * v[2968] - v[1231] * v[2969] - v[1230] * v[2970]
			- v[1229] * v[2971] - v[1228] * v[2972] - v[1227] * v[2973] - v[1226] * v[2974] - v[1225] * v[2975] - v[1224] * v[2976]
			- v[1223] * v[2977] - v[1222] * v[2978] - v[1221] * v[2979] - v[1220] * v[2980] - v[1219] * v[2981] - v[1218] * v[2982]
			- v[1217] * v[2983] - v[1216] * v[2984] - v[1215] * v[2985] - v[1214] * v[2986] - v[1213] * v[2987] - v[1212] * v[2988]
			- v[1211] * v[2989]) / 2e0;
		v[3430] = v[1258] * v[2966] + v[1257] * v[2967] + v[1256] * v[2968] + v[1255] * v[2969] + v[1254] * v[2970]
			+ v[1253] * v[2971] + v[1252] * v[2972] + v[1251] * v[2973] + v[1250] * v[2974] + v[1249] * v[2975] + v[1248] * v[2976]
			+ v[1247] * v[2977] + v[1246] * v[2978] + v[1245] * v[2979] + v[1244] * v[2980] + v[1243] * v[2981] + v[1242] * v[2982]
			+ v[1241] * v[2983] + v[1240] * v[2984] + v[1239] * v[2985] + v[1238] * v[2986] + v[1237] * v[2987] + v[1236] * v[2988]
			+ v[1235] * v[2989];
		v[5104] = v[239] * v[3430];
		v[5103] = v[240] * v[3430];
		v[3613] = (v[1282] * v[2917] + v[1281] * v[2918] + v[1280] * v[2919] + v[1279] * v[2920] + v[1278] * v[2921]
			+ v[1277] * v[2922] + v[1276] * v[2923] + v[1275] * v[2924] + v[1274] * v[2925] + v[1273] * v[2926] + v[1272] * v[2927]
			+ v[1271] * v[2928] + v[1270] * v[2929] + v[1269] * v[2930] + v[1268] * v[2931] + v[1267] * v[2932] + v[1266] * v[2933]
			+ v[1265] * v[2934] + v[1264] * v[2935] + v[1263] * v[2936] + v[1262] * v[2937] + v[1261] * v[2938] + v[1260] * v[2939]
			+ v[1259] * v[2940]) / 2e0;
		v[3432] = -(v[1306] * v[2917]) - v[1305] * v[2918] - v[1304] * v[2919] - v[1303] * v[2920] - v[1302] * v[2921]
			- v[1301] * v[2922] - v[1300] * v[2923] - v[1299] * v[2924] - v[1298] * v[2925] - v[1297] * v[2926] - v[1296] * v[2927]
			- v[1295] * v[2928] - v[1294] * v[2929] - v[1293] * v[2930] - v[1292] * v[2931] - v[1291] * v[2932] - v[1290] * v[2933]
			- v[1289] * v[2934] - v[1288] * v[2935] - v[1287] * v[2936] - v[1286] * v[2937] - v[1285] * v[2938] - v[1284] * v[2939]
			- v[1283] * v[2940];
		v[5084] = v[326] * v[3432];
		v[5083] = v[327] * v[3432];
		v[3747] = (-(v[1234] * v[2917]) - v[1233] * v[2918] - v[1232] * v[2919] - v[1231] * v[2920] - v[1230] * v[2921]
			- v[1229] * v[2922] - v[1228] * v[2923] - v[1227] * v[2924] - v[1226] * v[2925] - v[1225] * v[2926] - v[1224] * v[2927]
			- v[1223] * v[2928] - v[1222] * v[2929] - v[1221] * v[2930] - v[1220] * v[2931] - v[1219] * v[2932] - v[1218] * v[2933]
			- v[1217] * v[2934] - v[1216] * v[2935] - v[1215] * v[2936] - v[1214] * v[2937] - v[1213] * v[2938] - v[1212] * v[2939]
			- v[1211] * v[2940]) / 2e0;
		v[3434] = v[1258] * v[2917] + v[1257] * v[2918] + v[1256] * v[2919] + v[1255] * v[2920] + v[1254] * v[2921]
			+ v[1253] * v[2922] + v[1252] * v[2923] + v[1251] * v[2924] + v[1250] * v[2925] + v[1249] * v[2926] + v[1248] * v[2927]
			+ v[1247] * v[2928] + v[1246] * v[2929] + v[1245] * v[2930] + v[1244] * v[2931] + v[1243] * v[2932] + v[1242] * v[2933]
			+ v[1241] * v[2934] + v[1240] * v[2935] + v[1239] * v[2936] + v[1238] * v[2937] + v[1237] * v[2938] + v[1236] * v[2939]
			+ v[1235] * v[2940];
		v[5102] = v[239] * v[3434];
		v[5101] = v[240] * v[3434];
		v[3616] = (v[1282] * v[2868] + v[1281] * v[2869] + v[1280] * v[2870] + v[1279] * v[2871] + v[1278] * v[2872]
			+ v[1277] * v[2873] + v[1276] * v[2874] + v[1275] * v[2875] + v[1274] * v[2876] + v[1273] * v[2877] + v[1272] * v[2878]
			+ v[1271] * v[2879] + v[1270] * v[2880] + v[1269] * v[2881] + v[1268] * v[2882] + v[1267] * v[2883] + v[1266] * v[2884]
			+ v[1265] * v[2885] + v[1264] * v[2886] + v[1263] * v[2887] + v[1262] * v[2888] + v[1261] * v[2889] + v[1260] * v[2890]
			+ v[1259] * v[2891]) / 2e0;
		v[3436] = -(v[1306] * v[2868]) - v[1305] * v[2869] - v[1304] * v[2870] - v[1303] * v[2871] - v[1302] * v[2872]
			- v[1301] * v[2873] - v[1300] * v[2874] - v[1299] * v[2875] - v[1298] * v[2876] - v[1297] * v[2877] - v[1296] * v[2878]
			- v[1295] * v[2879] - v[1294] * v[2880] - v[1293] * v[2881] - v[1292] * v[2882] - v[1291] * v[2883] - v[1290] * v[2884]
			- v[1289] * v[2885] - v[1288] * v[2886] - v[1287] * v[2887] - v[1286] * v[2888] - v[1285] * v[2889] - v[1284] * v[2890]
			- v[1283] * v[2891];
		v[5082] = v[326] * v[3436];
		v[5081] = v[327] * v[3436];
		v[3750] = (-(v[1234] * v[2868]) - v[1233] * v[2869] - v[1232] * v[2870] - v[1231] * v[2871] - v[1230] * v[2872]
			- v[1229] * v[2873] - v[1228] * v[2874] - v[1227] * v[2875] - v[1226] * v[2876] - v[1225] * v[2877] - v[1224] * v[2878]
			- v[1223] * v[2879] - v[1222] * v[2880] - v[1221] * v[2881] - v[1220] * v[2882] - v[1219] * v[2883] - v[1218] * v[2884]
			- v[1217] * v[2885] - v[1216] * v[2886] - v[1215] * v[2887] - v[1214] * v[2888] - v[1213] * v[2889] - v[1212] * v[2890]
			- v[1211] * v[2891]) / 2e0;
		v[3438] = v[1258] * v[2868] + v[1257] * v[2869] + v[1256] * v[2870] + v[1255] * v[2871] + v[1254] * v[2872]
			+ v[1253] * v[2873] + v[1252] * v[2874] + v[1251] * v[2875] + v[1250] * v[2876] + v[1249] * v[2877] + v[1248] * v[2878]
			+ v[1247] * v[2879] + v[1246] * v[2880] + v[1245] * v[2881] + v[1244] * v[2882] + v[1243] * v[2883] + v[1242] * v[2884]
			+ v[1241] * v[2885] + v[1240] * v[2886] + v[1239] * v[2887] + v[1238] * v[2888] + v[1237] * v[2889] + v[1236] * v[2890]
			+ v[1235] * v[2891];
		v[5100] = v[239] * v[3438];
		v[5099] = v[240] * v[3438];
		v[3439] = -(v[1258] * v[2823]) - v[1234] * v[2824] + v[1306] * v[2825] + v[1282] * v[2829] + v[18] * (-
			(v[1658] * v[2836]) - v[1659] * v[2837] - v[1660] * v[2838]) + v[17] * (v[2835] * v[3065] + v[2834] * v[3167]
				+ v[2833] * v[3277]);
		v[3648] = v[3439] * v[4743];
		v[3440] = -(v[1257] * v[2823]) - v[1233] * v[2824] + v[1305] * v[2825] + v[1281] * v[2829] + v[18] * (-
			(v[1654] * v[2836]) - v[1655] * v[2837] - v[1656] * v[2838]) + v[17] * (v[2835] * v[3067] + v[2834] * v[3169]
				+ v[2833] * v[3279]);
		v[5053] = v[3440] * v[442];
		v[3652] = v[3440] * v[4742];
		v[3441] = -(v[1256] * v[2823]) - v[1232] * v[2824] + v[1304] * v[2825] + v[1280] * v[2829] + v[18] * (-
			(v[1650] * v[2836]) - v[1651] * v[2837] - v[1652] * v[2838]) + v[17] * (v[2835] * v[3069] + v[2834] * v[3171]
				+ v[2833] * v[3281]);
		v[5056] = v[3441] * v[442];
		v[3655] = v[3441] * v[4740];
		v[3442] = -(v[1252] * v[2823]) - v[1228] * v[2824] + v[1300] * v[2825] + v[1276] * v[2829] + v[18] * (-
			(v[1634] * v[2836]) - v[1635] * v[2837] - v[1636] * v[2838]) + v[17] * (v[2835] * v[3071] + v[2834] * v[3173]
				+ v[2833] * v[3283]);
		v[3678] = v[3442] * v[4739];
		v[3443] = -(v[1251] * v[2823]) - v[1227] * v[2824] + v[1299] * v[2825] + v[1275] * v[2829] + v[18] * (-
			(v[1630] * v[2836]) - v[1631] * v[2837] - v[1632] * v[2838]) + v[17] * (v[2835] * v[3073] + v[2834] * v[3175]
				+ v[2833] * v[3285]);
		v[5061] = v[3443] * v[416];
		v[3682] = v[3443] * v[4738];
		v[3444] = -(v[1250] * v[2823]) - v[1226] * v[2824] + v[1298] * v[2825] + v[1274] * v[2829] + v[18] * (-
			(v[1626] * v[2836]) - v[1627] * v[2837] - v[1628] * v[2838]) + v[17] * (v[2835] * v[3075] + v[2834] * v[3177]
				+ v[2833] * v[3287]);
		v[5064] = v[3444] * v[416];
		v[3685] = v[3444] * v[4736];
		v[3445] = -(v[1246] * v[2823]) - v[1222] * v[2824] + v[1294] * v[2825] + v[1270] * v[2829] + v[18] * (-
			(v[1610] * v[2836]) - v[1611] * v[2837] - v[1612] * v[2838]) + v[17] * (v[2835] * v[3077] + v[2834] * v[3179]
				+ v[2833] * v[3289]);
		v[3782] = v[3445] * v[4730];
		v[3446] = -(v[1245] * v[2823]) - v[1221] * v[2824] + v[1293] * v[2825] + v[1269] * v[2829] + v[18] * (-
			(v[1606] * v[2836]) - v[1607] * v[2837] - v[1608] * v[2838]) + v[17] * (v[2835] * v[3079] + v[2834] * v[3181]
				+ v[2833] * v[3291]);
		v[5069] = v[3446] * v[390];
		v[3786] = v[3446] * v[4729];
		v[3447] = -(v[1244] * v[2823]) - v[1220] * v[2824] + v[1292] * v[2825] + v[1268] * v[2829] + v[18] * (-
			(v[1602] * v[2836]) - v[1603] * v[2837] - v[1604] * v[2838]) + v[17] * (v[2835] * v[3081] + v[2834] * v[3183]
				+ v[2833] * v[3293]);
		v[5072] = v[3447] * v[390];
		v[3789] = v[3447] * v[4727];
		v[3448] = -(v[1240] * v[2823]) - v[1216] * v[2824] + v[1288] * v[2825] + v[1264] * v[2829] + v[18] * (-
			(v[1586] * v[2836]) - v[1587] * v[2837] - v[1588] * v[2838]) + v[17] * (v[2835] * v[3083] + v[2834] * v[3185]
				+ v[2833] * v[3295]);
		v[3812] = v[3448] * v[4726];
		v[3449] = -(v[1239] * v[2823]) - v[1215] * v[2824] + v[1287] * v[2825] + v[1263] * v[2829] + v[18] * (-
			(v[1582] * v[2836]) - v[1583] * v[2837] - v[1584] * v[2838]) + v[17] * (v[2835] * v[3085] + v[2834] * v[3187]
				+ v[2833] * v[3297]);
		v[5077] = v[3449] * v[364];
		v[3816] = v[3449] * v[4725];
		v[3450] = -(v[1238] * v[2823]) - v[1214] * v[2824] + v[1286] * v[2825] + v[1262] * v[2829] + v[18] * (-
			(v[1578] * v[2836]) - v[1579] * v[2837] - v[1580] * v[2838]) + v[17] * (v[2835] * v[3087] + v[2834] * v[3189]
				+ v[2833] * v[3299]);
		v[5080] = v[3450] * v[364];
		v[3819] = v[3450] * v[4723];
		v[2648] = v[2648] + 2e0 * v[3055] * v[482] + v[3242] * v[5046];
		v[2644] = v[2644] + 2e0 * v[3152] * v[481] + v[480] * v[5047];
		v[2640] = v[2640] + 2e0 * v[3250] * v[480] + v[481] * v[5047];
		v[3457] = v[175] * v[3324] + v[3323] * v[4801] + v[3317] * v[936] + v[3319] * v[944];
		v[3458] = v[175] * v[3321] + v[3319] * v[4800] + v[3317] * v[935] + v[3323] * v[950];
		v[3459] = v[175] * v[3322] + v[3317] * v[4799] + v[3319] * v[937] + v[3323] * v[945];
		v[3460] = v[194] * v[3336] + v[3335] * v[4798] + v[3329] * v[905] + v[3331] * v[913];
		v[3461] = v[194] * v[3333] + v[3331] * v[4797] + v[3329] * v[904] + v[3335] * v[919];
		v[3462] = v[194] * v[3334] + v[3329] * v[4796] + v[3331] * v[906] + v[3335] * v[914];
		v[3463] = v[262] * v[3348] + v[3347] * v[4789] + v[3341] * v[754] + v[3343] * v[762];
		v[3464] = v[262] * v[3345] + v[3343] * v[4788] + v[3341] * v[753] + v[3347] * v[768];
		v[3465] = v[262] * v[3346] + v[3341] * v[4787] + v[3343] * v[755] + v[3347] * v[763];
		v[3466] = v[281] * v[3360] + v[3359] * v[4786] + v[3353] * v[723] + v[3355] * v[731];
		v[3467] = v[281] * v[3357] + v[3355] * v[4785] + v[3353] * v[722] + v[3359] * v[737];
		v[3468] = v[281] * v[3358] + v[3353] * v[4784] + v[3355] * v[724] + v[3359] * v[732];
		v[3470] = v[1886] * v[2636] * v[477] + (v[1770] * v[2513] + v[1769] * v[2518] + v[1768] * v[2523] + v[2640] * v[465]
			+ v[2644] * v[466] + v[2648] * v[467]) * v[478] - v[1888] * v[3469] * v[5048] + (v[2833] * v[480] + v[2834] * v[481]
				+ v[2835] * v[482]) * v[5245];
		v[3472] = v[1770] * v[2637] + (v[1888] * v[2513] + v[3470] * v[467]) / v[1475] + v[2648] * v[479];
		v[3474] = v[1769] * v[2637] + (v[1888] * v[2518] + v[3470] * v[466]) / v[1475] + v[2644] * v[479];
		v[3476] = v[1768] * v[2637] + (v[1888] * v[2523] + v[3470] * v[465]) / v[1475] + v[2640] * v[479];
		v[3477] = v[3439] * v[451] + v[1737] * v[5049];
		v[3478] = v[3439] * v[445] + v[1737] * v[5050];
		v[3479] = v[1901] * v[3439] + (v[290] * v[3477]) / v[442] + v[1737] * (v[2485] / v[442] + v[2455] * v[5054]);
		v[3480] = v[1908] * v[3439] + (v[284] * v[3478]) / v[442] + v[1737] * (v[2480] / v[442] + v[2455] * v[5392]);
		v[3481] = (v[295] * v[3477] + v[293] * v[3478] + v[1897] * v[3439] * v[442] + v[1737] * (v[2490] + v[5003] * v[5393]))
			/ v[442];
		v[3482] = v[3440] * v[453] + v[1738] * v[5051];
		v[3483] = v[3440] * v[445] + v[1738] * v[5050];
		v[3486] = v[1896] * v[3440] + (v[296] * v[3482]) / v[442] + v[1738] * (v[2488] / v[442] + v[2455] * v[5052]);
		v[3487] = v[3477] + v[3482];
		v[3488] = -v[3479] + v[3486];
		v[3489] = v[3479] + v[3486];
		v[3491] = (v[1899] * v[2444] + v[284] * v[3483] + v[286] * v[3487] + v[3885] * v[5003] + v[1907] * v[5053] + v[1738] *
			(v[2476] + v[5003] * v[5055])) / v[442];
		v[3494] = (v[291] * v[3482] + v[288] * v[3483] + v[1900] * v[5053] + v[1738] * (v[2482] + v[5003] * v[5394])) / v[442];
		v[3495] = v[3441] * v[451] + v[1739] * v[5049];
		v[3496] = v[3441] * v[453] + v[1739] * v[5051];
		v[3497] = v[3483] + v[3495];
		v[3500] = (v[285] * v[3495] + v[286] * v[3496] + v[1905] * v[5056] + v[1739] * (v[2477] + v[5003] * v[5395])) / v[442];
		v[3501] = v[3478] + v[3496];
		v[3503] = (v[1910] * v[2447] + v[290] * v[3495] + v[291] * v[3501] + v[3884] * v[5003] + v[1912] * v[5056] + v[1739] *
			(v[2481] + v[5003] * v[5396])) / v[442];
		v[3504] = v[3491] + v[3503];
		v[3505] = -v[3491] + v[3503];
		v[3507] = (v[1911] * v[2449] + v[296] * v[3496] + v[295] * v[3497] + v[3883] * v[5003] + v[1916] * v[5056] + v[1739] *
			(v[2486] + v[5003] * v[5397])) / v[442];
		v[3508] = v[3480] + v[3507];
		v[3509] = -v[3480] + v[3507];
		v[3510] = v[3442] * v[425] + v[1740] * v[5057];
		v[3511] = v[3442] * v[419] + v[1740] * v[5058];
		v[3512] = v[1928] * v[3442] + (v[271] * v[3510]) / v[416] + v[1740] * (v[2470] / v[416] + v[2433] * v[5062]);
		v[3513] = v[1935] * v[3442] + (v[265] * v[3511]) / v[416] + v[1740] * (v[2465] / v[416] + v[2433] * v[5398]);
		v[3514] = (v[276] * v[3510] + v[274] * v[3511] + v[1924] * v[3442] * v[416] + v[1740] * (v[2475] + v[5002] * v[5399]))
			/ v[416];
		v[3515] = v[3443] * v[427] + v[1741] * v[5059];
		v[3516] = v[3443] * v[419] + v[1741] * v[5058];
		v[3519] = v[1923] * v[3443] + (v[277] * v[3515]) / v[416] + v[1741] * (v[2473] / v[416] + v[2433] * v[5060]);
		v[3520] = v[3510] + v[3515];
		v[3521] = -v[3512] + v[3519];
		v[3522] = v[3512] + v[3519];
		v[3524] = (v[1926] * v[2422] + v[265] * v[3516] + v[267] * v[3520] + v[3907] * v[5002] + v[1934] * v[5061] + v[1741] *
			(v[2461] + v[5002] * v[5063])) / v[416];
		v[3527] = (v[272] * v[3515] + v[269] * v[3516] + v[1927] * v[5061] + v[1741] * (v[2467] + v[5002] * v[5400])) / v[416];
		v[3528] = v[3444] * v[425] + v[1742] * v[5057];
		v[3529] = v[3444] * v[427] + v[1742] * v[5059];
		v[3530] = v[3516] + v[3528];
		v[3533] = (v[266] * v[3528] + v[267] * v[3529] + v[1932] * v[5064] + v[1742] * (v[2462] + v[5002] * v[5401])) / v[416];
		v[3534] = v[3511] + v[3529];
		v[3536] = (v[1937] * v[2425] + v[271] * v[3528] + v[272] * v[3534] + v[3906] * v[5002] + v[1939] * v[5064] + v[1742] *
			(v[2466] + v[5002] * v[5402])) / v[416];
		v[3537] = v[3524] + v[3536];
		v[3538] = -v[3524] + v[3536];
		v[3540] = (v[1938] * v[2427] + v[277] * v[3529] + v[276] * v[3530] + v[3905] * v[5002] + v[1943] * v[5064] + v[1742] *
			(v[2471] + v[5002] * v[5403])) / v[416];
		v[3541] = v[3513] + v[3540];
		v[3542] = -v[3513] + v[3540];
		v[3543] = v[3445] * v[399] + v[1743] * v[5065];
		v[3544] = v[3445] * v[393] + v[1743] * v[5066];
		v[3545] = v[1955] * v[3445] + (v[203] * v[3543]) / v[390] + v[1743] * (v[2381] / v[390] + v[2351] * v[5070]);
		v[3546] = v[1962] * v[3445] + (v[197] * v[3544]) / v[390] + v[1743] * (v[2376] / v[390] + v[2351] * v[5404]);
		v[3547] = (v[208] * v[3543] + v[206] * v[3544] + v[1951] * v[3445] * v[390] + v[1743] * (v[2386] + v[5001] * v[5405]))
			/ v[390];
		v[3548] = v[3446] * v[401] + v[1744] * v[5067];
		v[3549] = v[3446] * v[393] + v[1744] * v[5066];
		v[3552] = v[1950] * v[3446] + (v[209] * v[3548]) / v[390] + v[1744] * (v[2384] / v[390] + v[2351] * v[5068]);
		v[3553] = v[3543] + v[3548];
		v[3554] = -v[3545] + v[3552];
		v[3555] = v[3545] + v[3552];
		v[3557] = (v[1953] * v[2340] + v[197] * v[3549] + v[199] * v[3553] + v[3929] * v[5001] + v[1961] * v[5069] + v[1744] *
			(v[2372] + v[5001] * v[5071])) / v[390];
		v[3560] = (v[204] * v[3548] + v[201] * v[3549] + v[1954] * v[5069] + v[1744] * (v[2378] + v[5001] * v[5406])) / v[390];
		v[3561] = v[3447] * v[399] + v[1745] * v[5065];
		v[3562] = v[3447] * v[401] + v[1745] * v[5067];
		v[3563] = v[3549] + v[3561];
		v[3566] = (v[198] * v[3561] + v[199] * v[3562] + v[1959] * v[5072] + v[1745] * (v[2373] + v[5001] * v[5407])) / v[390];
		v[3567] = v[3544] + v[3562];
		v[3569] = (v[1964] * v[2343] + v[203] * v[3561] + v[204] * v[3567] + v[3928] * v[5001] + v[1966] * v[5072] + v[1745] *
			(v[2377] + v[5001] * v[5408])) / v[390];
		v[3570] = v[3557] + v[3569];
		v[3571] = -v[3557] + v[3569];
		v[3573] = (v[1965] * v[2345] + v[209] * v[3562] + v[208] * v[3563] + v[3927] * v[5001] + v[1970] * v[5072] + v[1745] *
			(v[2382] + v[5001] * v[5409])) / v[390];
		v[3574] = v[3546] + v[3573];
		v[3575] = -v[3546] + v[3573];
		v[3576] = v[3448] * v[373] + v[1746] * v[5073];
		v[3577] = v[3448] * v[367] + v[1746] * v[5074];
		v[3578] = v[1982] * v[3448] + (v[184] * v[3576]) / v[364] + v[1746] * (v[2366] / v[364] + v[2329] * v[5078]);
		v[3579] = v[1989] * v[3448] + (v[178] * v[3577]) / v[364] + v[1746] * (v[2361] / v[364] + v[2329] * v[5410]);
		v[3580] = (v[189] * v[3576] + v[187] * v[3577] + v[1978] * v[3448] * v[364] + v[1746] * (v[2371] + v[5000] * v[5411]))
			/ v[364];
		v[3581] = v[3449] * v[375] + v[1747] * v[5075];
		v[3582] = v[3449] * v[367] + v[1747] * v[5074];
		v[3585] = v[1977] * v[3449] + (v[190] * v[3581]) / v[364] + v[1747] * (v[2369] / v[364] + v[2329] * v[5076]);
		v[3586] = v[3576] + v[3581];
		v[3587] = -v[3578] + v[3585];
		v[3588] = v[3578] + v[3585];
		v[3590] = (v[1980] * v[2318] + v[178] * v[3582] + v[180] * v[3586] + v[3951] * v[5000] + v[1988] * v[5077] + v[1747] *
			(v[2357] + v[5000] * v[5079])) / v[364];
		v[3593] = (v[185] * v[3581] + v[182] * v[3582] + v[1981] * v[5077] + v[1747] * (v[2363] + v[5000] * v[5412])) / v[364];
		v[3594] = v[3450] * v[373] + v[1748] * v[5073];
		v[3595] = v[3450] * v[375] + v[1748] * v[5075];
		v[3596] = v[3582] + v[3594];
		v[3599] = (v[179] * v[3594] + v[180] * v[3595] + v[1986] * v[5080] + v[1748] * (v[2358] + v[5000] * v[5413])) / v[364];
		v[3600] = v[3577] + v[3595];
		v[3602] = (v[1991] * v[2321] + v[184] * v[3594] + v[185] * v[3600] + v[3950] * v[5000] + v[1993] * v[5080] + v[1748] *
			(v[2362] + v[5000] * v[5414])) / v[364];
		v[3603] = v[3590] + v[3602];
		v[3604] = -v[3590] + v[3602];
		v[3606] = (v[1992] * v[2323] + v[190] * v[3595] + v[189] * v[3596] + v[3949] * v[5000] + v[1997] * v[5080] + v[1748] *
			(v[2367] + v[5000] * v[5415])) / v[364];
		v[3607] = v[3579] + v[3606];
		v[3608] = -v[3579] + v[3606];
		v[3609] = -(v[326] * v[3476]) + v[3610];
		v[3611] = -(v[327] * v[3476]) - v[3610];
		v[3372] = v[3372] + v[3609];
		v[3371] = v[3371] + v[3611];
		v[3612] = -(v[326] * v[3474]) + v[3613];
		v[3614] = -(v[327] * v[3474]) - v[3613];
		v[3384] = v[3384] + v[3612];
		v[3383] = v[3383] + v[3614];
		v[3615] = -(v[326] * v[3472]) + v[3616];
		v[3617] = -(v[327] * v[3472]) - v[3616];
		v[3396] = v[3396] + v[3615];
		v[3618] = v[2001] * v[2498] + v[2003] * v[2500] + v[2004] * v[2502] + v[2006] * v[2504] + v[2007] * v[2506]
			+ v[2009] * v[2508] + v[326] * (-(v[149] * v[3463]) - v[152] * v[3464] - v[155] * v[3465]) + v[327] * (-(v[158] * v[3466]
				) - v[161] * v[3467] - v[164] * v[3468]) + v[301] * v[3609] + v[310] * v[3611] + v[304] * v[3612] + v[313] * v[3614]
			+ v[307] * v[3615] + v[316] * v[3617];
		v[3619] = v[2001] * v[2497] + v[2003] * v[2499] + v[2004] * v[2501] + v[2006] * v[2503] + v[2007] * v[2505]
			+ v[2009] * v[2507] + v[326] * (-(v[148] * v[3463]) - v[151] * v[3464] - v[154] * v[3465]) + v[327] * (-(v[157] * v[3466]
				) - v[160] * v[3467] - v[163] * v[3468]) + v[300] * v[3609] + v[309] * v[3611] + v[303] * v[3612] + v[312] * v[3614]
			+ v[306] * v[3615] + v[315] * v[3617];
		v[3395] = v[3395] + v[3617];
		v[3397] = v[3397] + v[339] * v[3617] + v[5081] * v[537];
		v[3398] = v[3398] + v[337] * v[3617] + v[5081] * v[536];
		v[3399] = v[3399] + v[339] * v[3615] + v[5082] * v[537];
		v[3400] = v[3400] + v[337] * v[3615] + v[5082] * v[536];
		v[3385] = v[3385] + v[339] * v[3614] + v[5083] * v[537];
		v[3386] = v[3386] + v[337] * v[3614] + v[5083] * v[536];
		v[3387] = v[3387] + v[339] * v[3612] + v[5084] * v[537];
		v[3388] = v[3388] + v[337] * v[3612] + v[5084] * v[536];
		v[3620] = (v[1762] * v[2497] + v[1731] * v[2501] + v[1700] * v[2505]) * v[326] + (v[1762] * v[2499] + v[1731] * v[2503]
			+ v[1700] * v[2507]) * v[327] + v[3428] * v[4032] + v[3432] * v[4033] + v[3436] * v[4034];
		v[5165] = v[139] * v[3620];
		v[3621] = (v[1762] * v[2498] + v[1731] * v[2502] + v[1700] * v[2506]) * v[326] + (v[1762] * v[2500] + v[1731] * v[2504]
			+ v[1700] * v[2508]) * v[327] + v[3428] * v[4028] + v[3432] * v[4029] + v[3436] * v[4030];
		v[5166] = v[140] * v[3621];
		v[3373] = v[3373] + v[339] * v[3611] + v[5085] * v[537];
		v[3374] = v[3374] + v[337] * v[3611] + v[5085] * v[536];
		v[3375] = v[3375] + v[339] * v[3609] + v[5086] * v[537];
		v[3376] = v[3376] + v[337] * v[3609] + v[5086] * v[536];
		v[3622] = v[3621] * v[4748] + v[3620] * v[5087];
		v[5167] = v[328] * ((v[139] * (v[1565] * v[3619] - v[1462] * v[3620]) + v[140] * (v[1462] * v[3618] + v[1565] * v[3621])
			) * v[5088] + v[3622] * v[5265]);
		v[3632] = (-(v[1889] * v[2509]) + v[1889] * v[2511] - v[1891] * v[2514] + v[1891] * v[2516] - v[1892] * v[2519]
			+ v[1892] * v[2521] + v[3212] * v[355] + v[3112] * v[356] - v[3212] * v[358] - v[3112] * v[359] + v[3463] * v[3626]
			+ v[3464] * v[3627] + v[3465] * v[3628] - v[3466] * v[3629] - v[3467] * v[3630] - v[3468] * v[3631] - v[4526] * v[5028]
			+ v[4527] * v[5028] + v[357] * v[5089] - v[360] * v[5089] + v[3248] * v[5090] - v[3249] * v[5090] + v[3476] * v[520]
			- v[3476] * v[521] + v[3474] * v[523] - v[3474] * v[524] + v[3472] * v[526] - v[3472] * v[527] + (v[1762] * (-v[2497]
				+ v[2499]) + v[1731] * (-v[2501] + v[2503]) + v[1700] * (-v[2505] + v[2507])) * v[536] + (v[1762] * (-v[2498] + v[2500]
					) + v[1731] * (-v[2502] + v[2504]) + v[1700] * (-v[2506] + v[2508])) * v[537] + v[3436] * v[538] - v[3436] * v[539]
			+ v[3432] * v[540] - v[3432] * v[541] + v[3428] * v[542] - v[3428] * v[543]) / 2e0;
		v[3699] = v[1918] * v[2439] - v[1914] * v[2440] + 2e0 * v[1909] * v[2441] - v[281] * v[3488] - v[298] * v[3504]
			+ v[297] * v[3508] - v[1903] * v[4990] + v[3500] * v[4993];
		v[3700] = -(v[1904] * v[2439]) + 2e0 * v[1902] * v[2440] - v[1914] * v[2441] - v[297] * v[3489] - v[299] * v[3504]
			- v[281] * v[3509] - v[1919] * v[4990] + v[3494] * v[4992];
		v[3701] = 2e0 * v[1895] * v[2439] - v[1904] * v[2440] + v[1918] * v[2441] - v[298] * v[3489] - v[281] * v[3505]
			+ v[299] * v[3508] - v[1915] * v[4990] + v[3481] * v[4991];
		v[3702] = v[164] * v[3397] + v[163] * v[3398] + v[2227] * v[3482] + v[3648] * v[455] + v[3496] * v[4767] + v[453] *
			(v[5222] + v[2455] * v[5416]);
		v[3703] = v[2022] * v[3711] + v[3702] * v[4741] + v[3353] * v[5274];
		v[3704] = v[161] * v[3397] + v[160] * v[3398] + v[3477] * v[3998] + (v[1911] * v[2608]) / v[442] + v[3652] * v[452]
			+ v[2625] * v[4742] + v[3497] * v[4767] + v[2455] * v[5417];
		v[3705] = v[1348] * v[3353] + v[281] * v[3704] + v[2023] * v[4990];
		v[3706] = v[158] * v[3397] + v[157] * v[3398] + v[3478] * v[3998] + v[3655] * v[440] + v[445] * (v[5215]
			+ v[2455] * v[5418]);
		v[3707] = v[1349] * v[3353] + v[281] * v[3706] + v[2024] * v[4990];
		v[3708] = v[164] * v[3385] + v[163] * v[3386] + v[3482] * v[3996] + (v[1910] * v[2612]) / v[442] + v[2633] * v[4743]
			+ v[3501] * v[4765] + v[3648] * v[4768] + v[2455] * v[5419];
		v[3709] = v[1347] * v[3355] + v[281] * v[3708] + v[2031] * v[4990];
		v[3710] = v[161] * v[3385] + v[160] * v[3386] + v[2229] * v[3477] + v[3652] * v[446] + v[3495] * v[4765] + v[451] *
			(v[5219] + v[2455] * v[5420]);
		v[3713] = v[158] * v[3385] + v[157] * v[3386] + v[3483] * v[3996] + v[3655] * v[441] + v[445] * (v[5214]
			+ v[2455] * v[5421]);
		v[3714] = v[1349] * v[3355] + v[281] * v[3713] + v[2033] * v[4990];
		v[3715] = v[164] * v[3373] + v[163] * v[3374] + v[2228] * v[3487] + v[3496] * v[3994] + (v[1899] * v[2613]) / v[442]
			+ v[2629] * v[4743] + v[3648] * v[4766] + v[2455] * v[5422];
		v[3716] = v[1347] * v[3359] + v[281] * v[3715] + v[2040] * v[4990];
		v[3717] = v[161] * v[3373] + v[160] * v[3374] + v[3495] * v[3994] + v[3652] * v[4764] + v[451] * (v[2455] * v[5095]
			+ v[5218]);
		v[3718] = v[1348] * v[3359] + v[281] * v[3717] + v[2041] * v[4990];
		v[3719] = v[158] * v[3373] + v[157] * v[3374] + v[2230] * v[3478] + v[2228] * v[3483] + v[3655] * v[443] + v[445] *
			(v[5216] + v[2455] * v[5423]);
		v[3720] = v[2042] * v[3711] + v[3719] * v[4741] + v[3359] * v[5276];
		v[3721] = v[1945] * v[2417] - v[1941] * v[2418] + 2e0 * v[1936] * v[2419] - v[262] * v[3521] - v[279] * v[3537]
			+ v[278] * v[3541] - v[1930] * v[4986] + v[3533] * v[4989];
		v[3722] = -(v[1931] * v[2417]) + 2e0 * v[1929] * v[2418] - v[1941] * v[2419] - v[278] * v[3522] - v[280] * v[3537]
			- v[262] * v[3542] - v[1946] * v[4986] + v[3527] * v[4988];
		v[3723] = 2e0 * v[1922] * v[2417] - v[1931] * v[2418] + v[1945] * v[2419] - v[279] * v[3522] - v[262] * v[3538]
			+ v[280] * v[3541] - v[1942] * v[4986] + v[3514] * v[4987];
		v[3724] = v[155] * v[3399] + v[154] * v[3400] + v[2221] * v[3515] + v[3678] * v[429] + v[3529] * v[4762] + v[427] *
			(v[5209] + v[2433] * v[5424]);
		v[3725] = v[2049] * v[3733] + v[3724] * v[4737] + v[3341] * v[5277];
		v[3726] = v[152] * v[3399] + v[151] * v[3400] + v[3510] * v[3992] + (v[1938] * v[2580]) / v[416] + v[3682] * v[426]
			+ v[2597] * v[4738] + v[3530] * v[4762] + v[2433] * v[5425];
		v[3727] = v[1351] * v[3341] + v[262] * v[3726] + v[2050] * v[4986];
		v[3728] = v[149] * v[3399] + v[148] * v[3400] + v[3511] * v[3992] + v[3685] * v[414] + v[419] * (v[5202]
			+ v[2433] * v[5426]);
		v[3729] = v[1352] * v[3341] + v[262] * v[3728] + v[2051] * v[4986];
		v[3730] = v[155] * v[3387] + v[154] * v[3388] + v[3515] * v[3990] + (v[1937] * v[2584]) / v[416] + v[2605] * v[4739]
			+ v[3534] * v[4760] + v[3678] * v[4763] + v[2433] * v[5427];
		v[3731] = v[1350] * v[3343] + v[262] * v[3730] + v[2058] * v[4986];
		v[3732] = v[152] * v[3387] + v[151] * v[3388] + v[2223] * v[3510] + v[3682] * v[420] + v[3528] * v[4760] + v[425] *
			(v[5206] + v[2433] * v[5428]);
		v[3735] = v[149] * v[3387] + v[148] * v[3388] + v[3516] * v[3990] + v[3685] * v[415] + v[419] * (v[5201]
			+ v[2433] * v[5429]);
		v[3736] = v[1352] * v[3343] + v[262] * v[3735] + v[2060] * v[4986];
		v[3737] = v[155] * v[3375] + v[154] * v[3376] + v[2222] * v[3520] + v[3529] * v[3988] + (v[1926] * v[2585]) / v[416]
			+ v[2601] * v[4739] + v[3678] * v[4761] + v[2433] * v[5430];
		v[3738] = v[1350] * v[3347] + v[262] * v[3737] + v[2067] * v[4986];
		v[3739] = v[152] * v[3375] + v[151] * v[3376] + v[3528] * v[3988] + v[3682] * v[4759] + v[425] * (v[2433] * v[5098]
			+ v[5205]);
		v[3740] = v[1351] * v[3347] + v[262] * v[3739] + v[2068] * v[4986];
		v[3741] = v[149] * v[3375] + v[148] * v[3376] + v[2224] * v[3511] + v[2222] * v[3516] + v[3685] * v[417] + v[419] *
			(v[5203] + v[2433] * v[5431]);
		v[3742] = v[2069] * v[3733] + v[3741] * v[4737] + v[3347] * v[5279];
		v[3743] = v[239] * v[3476] + v[3744];
		v[3745] = v[240] * v[3476] - v[3744];
		v[3378] = v[3378] + v[3743];
		v[3377] = v[3377] + v[3745];
		v[3746] = v[239] * v[3474] + v[3747];
		v[3748] = v[240] * v[3474] - v[3747];
		v[3390] = v[3390] + v[3746];
		v[3389] = v[3389] + v[3748];
		v[3749] = v[239] * v[3472] + v[3750];
		v[3751] = v[240] * v[3472] - v[3750];
		v[3402] = v[3402] + v[3749];
		v[3752] = v[2098] * v[2394] + v[2100] * v[2396] + v[2101] * v[2398] + v[2103] * v[2400] + v[2104] * v[2402]
			+ v[2106] * v[2404] + v[214] * v[3743] + v[223] * v[3745] + v[217] * v[3746] + v[226] * v[3748] + v[220] * v[3749]
			+ v[229] * v[3751] + v[239] * (v[3457] * v[86] + v[3458] * v[89] + v[3459] * v[92]) + v[240] * (v[101] * v[3462]
				+ v[3460] * v[95] + v[3461] * v[98]);
		v[3753] = v[2098] * v[2393] + v[2100] * v[2395] + v[2101] * v[2397] + v[2103] * v[2399] + v[2104] * v[2401]
			+ v[2106] * v[2403] + v[213] * v[3743] + v[222] * v[3745] + v[216] * v[3746] + v[225] * v[3748] + v[219] * v[3749]
			+ v[228] * v[3751] + v[239] * (v[3457] * v[85] + v[3458] * v[88] + v[3459] * v[91]) + v[240] * (v[100] * v[3462]
				+ v[3460] * v[94] + v[3461] * v[97]);
		v[3401] = v[3401] + v[3751];
		v[3403] = v[3403] + v[252] * v[3751] + v[509] * v[5099];
		v[3404] = v[3404] + v[250] * v[3751] + v[508] * v[5099];
		v[3405] = v[3405] + v[252] * v[3749] + v[509] * v[5100];
		v[3406] = v[3406] + v[250] * v[3749] + v[508] * v[5100];
		v[3391] = v[3391] + v[252] * v[3748] + v[509] * v[5101];
		v[3392] = v[3392] + v[250] * v[3748] + v[508] * v[5101];
		v[3393] = v[3393] + v[252] * v[3746] + v[509] * v[5102];
		v[3394] = v[3394] + v[250] * v[3746] + v[508] * v[5102];
		v[3754] = v[239] * (v[1764] * v[2393] + v[1735] * v[2397] + v[1702] * v[2401]) + v[240] * (v[1764] * v[2395]
			+ v[1735] * v[2399] + v[1702] * v[2403]) + v[3430] * v[4046] + v[3434] * v[4047] + v[3438] * v[4048];
		v[5168] = v[3754] * v[76];
		v[3755] = v[239] * (v[1764] * v[2394] + v[1735] * v[2398] + v[1702] * v[2402]) + v[240] * (v[1764] * v[2396]
			+ v[1735] * v[2400] + v[1702] * v[2404]) + v[3430] * v[4042] + v[3434] * v[4043] + v[3438] * v[4044];
		v[5169] = v[3755] * v[77];
		v[3379] = v[3379] + v[252] * v[3745] + v[509] * v[5103];
		v[3380] = v[3380] + v[250] * v[3745] + v[508] * v[5103];
		v[3381] = v[3381] + v[252] * v[3743] + v[509] * v[5104];
		v[3382] = v[3382] + v[250] * v[3743] + v[508] * v[5104];
		v[3756] = v[3755] * v[4735] + v[3754] * v[5105];
		v[5170] = v[241] * (v[3756] * v[5286] + v[5106] * ((v[1564] * v[3753] - v[1457] * v[3754]) * v[76] + (v[1457] * v[3752]
			+ v[1564] * v[3755]) * v[77]));
		v[3766] = (v[1889] * v[2405] - v[1889] * v[2407] + v[1891] * v[2409] - v[1891] * v[2411] + v[1892] * v[2413]
			- v[1892] * v[2415] - v[3212] * v[349] - v[3112] * v[350] + v[3212] * v[352] + v[3112] * v[353] - v[3457] * v[3760]
			- v[3458] * v[3761] - v[3459] * v[3762] + v[3460] * v[3763] + v[3461] * v[3764] + v[3462] * v[3765] - v[3476] * v[492]
			+ v[3476] * v[493] - v[3474] * v[495] + v[3474] * v[496] - v[3472] * v[498] + v[3472] * v[499] - v[4607] * v[5028]
			+ v[4608] * v[5028] + (v[1764] * (-v[2393] + v[2395]) + v[1735] * (-v[2397] + v[2399]) + v[1702] * (-v[2401] + v[2403])
				) * v[508] - v[351] * v[5089] + v[354] * v[5089] + (v[1764] * (-v[2394] + v[2396]) + v[1735] * (-v[2398] + v[2400])
					+ v[1702] * (-v[2402] + v[2404])) * v[509] - v[3246] * v[5090] + v[3247] * v[5090] + v[3438] * v[510] - v[3438] * v[511]
			+ v[3434] * v[512] - v[3434] * v[513] + v[3430] * v[514] - v[3430] * v[515]) / 2e0;
		v[3833] = v[1972] * v[2335] - v[1968] * v[2336] + 2e0 * v[1963] * v[2337] - v[194] * v[3554] - v[211] * v[3570]
			+ v[210] * v[3574] - v[1957] * v[4976] + v[3566] * v[4979];
		v[3834] = -(v[1958] * v[2335]) + 2e0 * v[1956] * v[2336] - v[1968] * v[2337] - v[210] * v[3555] - v[212] * v[3570]
			- v[194] * v[3575] - v[1973] * v[4976] + v[3560] * v[4978];
		v[3835] = 2e0 * v[1949] * v[2335] - v[1958] * v[2336] + v[1972] * v[2337] - v[211] * v[3555] - v[194] * v[3571]
			+ v[212] * v[3574] - v[1969] * v[4976] + v[3547] * v[4977];
		v[3836] = v[101] * v[3403] + v[100] * v[3404] + v[2215] * v[3548] + v[3782] * v[403] + v[3562] * v[4757] + v[401] *
			(v[5196] + v[2351] * v[5432]);
		v[3837] = v[2119] * v[3845] + v[3836] * v[4728] + v[3329] * v[5293];
		v[3838] = (v[1965] * v[2552]) / v[390] + v[3543] * v[3986] + v[3786] * v[400] + v[2569] * v[4729] + v[3563] * v[4757]
			+ v[2351] * v[5433] + v[3404] * v[97] + v[3403] * v[98];
		v[3839] = v[1354] * v[3329] + v[194] * v[3838] + v[2120] * v[4976];
		v[3840] = v[3789] * v[388] + v[3544] * v[3986] + v[393] * (v[5189] + v[2351] * v[5434]) + v[3404] * v[94]
			+ v[3403] * v[95];
		v[3841] = v[1355] * v[3329] + v[194] * v[3840] + v[2121] * v[4976];
		v[3842] = v[101] * v[3391] + v[100] * v[3392] + (v[1964] * v[2556]) / v[390] + v[3548] * v[3984] + v[2577] * v[4730]
			+ v[3567] * v[4755] + v[3782] * v[4758] + v[2351] * v[5435];
		v[3843] = v[1353] * v[3331] + v[194] * v[3842] + v[2128] * v[4976];
		v[3844] = v[2217] * v[3543] + v[3786] * v[394] + v[3561] * v[4755] + v[399] * (v[5193] + v[2351] * v[5436])
			+ v[3392] * v[97] + v[3391] * v[98];
		v[3847] = v[3789] * v[389] + v[3549] * v[3984] + v[393] * (v[5188] + v[2351] * v[5437]) + v[3392] * v[94]
			+ v[3391] * v[95];
		v[3848] = v[1355] * v[3331] + v[194] * v[3847] + v[2130] * v[4976];
		v[3849] = v[101] * v[3379] + v[100] * v[3380] + v[2216] * v[3553] + (v[1953] * v[2557]) / v[390] + v[3562] * v[3982]
			+ v[2573] * v[4730] + v[3782] * v[4756] + v[2351] * v[5438];
		v[3850] = v[1353] * v[3335] + v[194] * v[3849] + v[2137] * v[4976];
		v[3851] = v[3561] * v[3982] + v[3786] * v[4754] + v[399] * (v[2351] * v[5109] + v[5192]) + v[3380] * v[97]
			+ v[3379] * v[98];
		v[3852] = v[1354] * v[3335] + v[194] * v[3851] + v[2138] * v[4976];
		v[3853] = v[2218] * v[3544] + v[2216] * v[3549] + v[3789] * v[391] + v[393] * (v[5190] + v[2351] * v[5439])
			+ v[3380] * v[94] + v[3379] * v[95];
		v[3854] = v[2139] * v[3845] + v[3853] * v[4728] + v[3335] * v[5295];
		v[3855] = v[1999] * v[2313] - v[1995] * v[2314] + 2e0 * v[1990] * v[2315] - v[175] * v[3587] - v[192] * v[3603]
			+ v[191] * v[3607] - v[1984] * v[4972] + v[3599] * v[4975];
		v[3856] = -(v[1985] * v[2313]) + 2e0 * v[1983] * v[2314] - v[1995] * v[2315] - v[191] * v[3588] - v[193] * v[3603]
			- v[175] * v[3608] - v[2000] * v[4972] + v[3593] * v[4974];
		v[3857] = 2e0 * v[1976] * v[2313] - v[1985] * v[2314] + v[1999] * v[2315] - v[192] * v[3588] - v[175] * v[3604]
			+ v[193] * v[3607] - v[1996] * v[4972] + v[3580] * v[4973];
		v[3858] = v[2209] * v[3581] + v[377] * v[3812] + v[3595] * v[4752] + v[375] * (v[5183] + v[2329] * v[5440])
			+ v[3406] * v[91] + v[3405] * v[92];
		v[3859] = v[2146] * v[3867] + v[3858] * v[4724] + v[3317] * v[5296];
		v[3860] = (v[1992] * v[2524]) / v[364] + v[374] * v[3816] + v[3576] * v[3980] + v[2541] * v[4725] + v[3596] * v[4752]
			+ v[2329] * v[5441] + v[3406] * v[88] + v[3405] * v[89];
		v[3861] = v[1357] * v[3317] + v[175] * v[3860] + v[2147] * v[4972];
		v[3862] = v[362] * v[3819] + v[3577] * v[3980] + v[367] * (v[5176] + v[2329] * v[5442]) + v[3406] * v[85]
			+ v[3405] * v[86];
		v[3863] = v[1358] * v[3317] + v[175] * v[3862] + v[2148] * v[4972];
		v[3864] = (v[1991] * v[2528]) / v[364] + v[3581] * v[3978] + v[2549] * v[4726] + v[3600] * v[4750] + v[3812] * v[4753]
			+ v[2329] * v[5443] + v[3394] * v[91] + v[3393] * v[92];
		v[3865] = v[1356] * v[3319] + v[175] * v[3864] + v[2155] * v[4972];
		v[3866] = v[2211] * v[3576] + v[368] * v[3816] + v[3594] * v[4750] + v[373] * (v[5180] + v[2329] * v[5444])
			+ v[3394] * v[88] + v[3393] * v[89];
		v[3869] = v[363] * v[3819] + v[3582] * v[3978] + v[367] * (v[5175] + v[2329] * v[5445]) + v[3394] * v[85]
			+ v[3393] * v[86];
		v[3870] = v[1358] * v[3319] + v[175] * v[3869] + v[2157] * v[4972];
		v[3871] = v[2210] * v[3586] + (v[1980] * v[2529]) / v[364] + v[3595] * v[3976] + v[2545] * v[4726] + v[3812] * v[4751]
			+ v[2329] * v[5446] + v[3382] * v[91] + v[3381] * v[92];
		v[3872] = v[1356] * v[3323] + v[175] * v[3871] + v[2164] * v[4972];
		v[3873] = v[3594] * v[3976] + v[3816] * v[4749] + v[373] * (v[2329] * v[5112] + v[5179]) + v[3382] * v[88]
			+ v[3381] * v[89];
		v[3874] = v[1357] * v[3323] + v[175] * v[3873] + v[2165] * v[4972];
		v[3875] = v[2212] * v[3577] + v[2210] * v[3582] + v[365] * v[3819] + v[367] * (v[5177] + v[2329] * v[5447])
			+ v[3382] * v[85] + v[3381] * v[86];
		v[3876] = v[2166] * v[3867] + v[3875] * v[4724] + v[3323] * v[5298];
		v[3895] = (v[2267] * v[3877] + v[1364] * v[5123] + v[1363] * v[5124] + v[1362] * v[5125]) * v[5299] + v[749] *
			(v[2041] * v[2307] + v[2033] * v[2308] + v[2031] * v[2309] + v[2023] * v[2310] + v[2040] * v[2311] + v[2024] * v[2312]
				- v[1915] * v[2439] - v[1919] * v[2440] - v[1903] * v[2441] + v[1348] * v[3357] + v[1347] * v[3358] + v[1349] * v[3360]
				- v[299] * v[3488] - v[297] * v[3505] - v[298] * v[3509] + v[3878] * v[3879] + v[3702] * v[4784] + v[3710] * v[4785]
				+ v[3719] * v[4786] + v[5113] * (v[3481] + v[2625] * v[3484] + v[2622] * v[3485] + v[2616] * v[3487] + v[2631] * v[3490]
					+ v[2635] * v[3492] + v[2633] * v[3493] + v[3494] + v[2610] * v[3497] + v[2618] * v[3498] + v[2629] * v[3499] + v[3500]
					+ v[2614] * v[3501] + v[2621] * v[3502] + v[2627] * v[3506] + v[2490] * v[3640] + v[2488] * v[3641] + v[2486] * v[3643]
					+ v[2485] * v[3650] + v[2482] * v[3651] + v[2481] * v[3654] + v[2477] * v[3661] + v[2480] * v[3663] + v[2476] * v[3664]
					+ v[2444] * v[3880] + v[2447] * v[3881] + v[2449] * v[3882] + v[2608] * v[3883] + v[2612] * v[3884] + v[2613] * v[3885]
					+ v[3477] * v[3889] + v[3478] * v[3890] + v[3482] * v[3891] + v[3483] * v[3892] + v[3495] * v[3893] + v[3496] * v[3894]
					+ v[3439] * v[5120] + v[3440] * v[5121] + v[3441] * v[5122] - 2e0 * v[2455] * v[5449] * v[5450]) + v[3704] * v[722]
				+ v[3706] * v[723] + v[3708] * v[724] + v[3713] * v[731] + v[3715] * v[732] + v[3717] * v[737] + (v[2042] * v[2268]
					+ v[2032] * v[2269] + v[2022] * v[2270] - v[8629 + i1663]) / 2e0 - v[3699] * v[964] - v[3701] * v[965] + v[3700] * v[966]
				);
		v[5210] = -(v[2032] * v[3711]) + v[3895] - v[3710] * v[4741] - v[3355] * v[5275];
		v[3896] = v[3707] + v[3716];
		v[3897] = v[3705] + v[3709];
		v[3898] = v[3714] + v[3718];
		v[3917] = (v[2263] * v[3899] + v[1370] * v[5136] + v[1369] * v[5137] + v[1368] * v[5138]) * v[5303] + v[780] *
			(v[2068] * v[2301] + v[2060] * v[2302] + v[2058] * v[2303] + v[2050] * v[2304] + v[2067] * v[2305] + v[2051] * v[2306]
				- v[1942] * v[2417] - v[1946] * v[2418] - v[1930] * v[2419] + v[1351] * v[3345] + v[1350] * v[3346] + v[1352] * v[3348]
				- v[280] * v[3521] - v[278] * v[3538] - v[279] * v[3542] + v[3900] * v[3901] + v[3724] * v[4787] + v[3732] * v[4788]
				+ v[3741] * v[4789] + v[5126] * (v[3514] + v[2597] * v[3517] + v[2594] * v[3518] + v[2588] * v[3520] + v[2603] * v[3523]
					+ v[2607] * v[3525] + v[2605] * v[3526] + v[3527] + v[2582] * v[3530] + v[2590] * v[3531] + v[2601] * v[3532] + v[3533]
					+ v[2586] * v[3534] + v[2593] * v[3535] + v[2599] * v[3539] + v[2475] * v[3670] + v[2473] * v[3671] + v[2471] * v[3673]
					+ v[2470] * v[3680] + v[2467] * v[3681] + v[2466] * v[3684] + v[2462] * v[3691] + v[2465] * v[3693] + v[2461] * v[3694]
					+ v[2422] * v[3902] + v[2425] * v[3903] + v[2427] * v[3904] + v[2580] * v[3905] + v[2584] * v[3906] + v[2585] * v[3907]
					+ v[3510] * v[3911] + v[3511] * v[3912] + v[3515] * v[3913] + v[3516] * v[3914] + v[3528] * v[3915] + v[3529] * v[3916]
					+ v[3442] * v[5133] + v[3443] * v[5134] + v[3444] * v[5135] - 2e0 * v[2433] * v[5452] * v[5453]) + v[3726] * v[753]
				+ v[3728] * v[754] + v[3730] * v[755] + v[3735] * v[762] + v[3737] * v[763] + v[3739] * v[768] + (v[2069] * v[2264]
					+ v[2059] * v[2265] + v[2049] * v[2266] - v[8605 + i1663]) / 2e0 - v[3721] * v[970] - v[3723] * v[971] + v[3722] * v[972]
				);
		v[5197] = -(v[2059] * v[3733]) + v[3917] - v[3732] * v[4737] - v[3343] * v[5278];
		v[3918] = v[3729] + v[3738];
		v[3919] = v[3727] + v[3731];
		v[3920] = v[3736] + v[3740];
		v[3939] = (v[2253] * v[3921] + v[1376] * v[5149] + v[1375] * v[5150] + v[1374] * v[5151]) * v[5307] + v[931] *
			(v[2138] * v[2295] + v[2130] * v[2296] + v[2128] * v[2297] + v[2120] * v[2298] + v[2137] * v[2299] + v[2121] * v[2300]
				- v[1969] * v[2335] - v[1973] * v[2336] - v[1957] * v[2337] + v[1354] * v[3333] + v[1353] * v[3334] + v[1355] * v[3336]
				- v[212] * v[3554] - v[210] * v[3571] - v[211] * v[3575] + v[3922] * v[3923] + v[3836] * v[4796] + v[3844] * v[4797]
				+ v[3853] * v[4798] + v[5139] * (v[3547] + v[2569] * v[3550] + v[2566] * v[3551] + v[2560] * v[3553] + v[2575] * v[3556]
					+ v[2579] * v[3558] + v[2577] * v[3559] + v[3560] + v[2554] * v[3563] + v[2562] * v[3564] + v[2573] * v[3565] + v[3566]
					+ v[2558] * v[3567] + v[2565] * v[3568] + v[2571] * v[3572] + v[2386] * v[3774] + v[2384] * v[3775] + v[2382] * v[3777]
					+ v[2381] * v[3784] + v[2378] * v[3785] + v[2377] * v[3788] + v[2373] * v[3795] + v[2376] * v[3797] + v[2372] * v[3798]
					+ v[2340] * v[3924] + v[2343] * v[3925] + v[2345] * v[3926] + v[2552] * v[3927] + v[2556] * v[3928] + v[2557] * v[3929]
					+ v[3543] * v[3933] + v[3544] * v[3934] + v[3548] * v[3935] + v[3549] * v[3936] + v[3561] * v[3937] + v[3562] * v[3938]
					+ v[3445] * v[5146] + v[3446] * v[5147] + v[3447] * v[5148] - 2e0 * v[2351] * v[5455] * v[5456]) + (v[2139] * v[2254]
						+ v[2129] * v[2255] + v[2119] * v[2256] - v[8581 + i1663]) / 2e0 + v[3838] * v[904] + v[3840] * v[905] + v[3842] * v[906]
				+ v[3847] * v[913] + v[3849] * v[914] + v[3851] * v[919] - v[3833] * v[976] - v[3835] * v[977] + v[3834] * v[978]);
		v[5184] = -(v[2129] * v[3845]) + v[3939] - v[3844] * v[4728] - v[3331] * v[5294];
		v[3940] = v[3841] + v[3850];
		v[3941] = v[3839] + v[3843];
		v[3942] = v[3848] + v[3852];
		v[3961] = (v[2249] * v[3943] + v[1382] * v[5162] + v[1381] * v[5163] + v[1380] * v[5164]) * v[5311] + v[962] *
			(v[2165] * v[2289] + v[2157] * v[2290] + v[2155] * v[2291] + v[2147] * v[2292] + v[2164] * v[2293] + v[2148] * v[2294]
				- v[1996] * v[2313] - v[2000] * v[2314] - v[1984] * v[2315] + v[1357] * v[3321] + v[1356] * v[3322] + v[1358] * v[3324]
				- v[193] * v[3587] - v[191] * v[3604] - v[192] * v[3608] + v[3944] * v[3945] + v[3858] * v[4799] + v[3866] * v[4800]
				+ v[3875] * v[4801] + v[5152] * (v[3580] + v[2541] * v[3583] + v[2538] * v[3584] + v[2532] * v[3586] + v[2547] * v[3589]
					+ v[2551] * v[3591] + v[2549] * v[3592] + v[3593] + v[2526] * v[3596] + v[2534] * v[3597] + v[2545] * v[3598] + v[3599]
					+ v[2530] * v[3600] + v[2537] * v[3601] + v[2543] * v[3605] + v[2371] * v[3804] + v[2369] * v[3805] + v[2367] * v[3807]
					+ v[2366] * v[3814] + v[2363] * v[3815] + v[2362] * v[3818] + v[2358] * v[3825] + v[2361] * v[3827] + v[2357] * v[3828]
					+ v[2318] * v[3946] + v[2321] * v[3947] + v[2323] * v[3948] + v[2524] * v[3949] + v[2528] * v[3950] + v[2529] * v[3951]
					+ v[3576] * v[3955] + v[3577] * v[3956] + v[3581] * v[3957] + v[3582] * v[3958] + v[3594] * v[3959] + v[3595] * v[3960]
					+ v[3448] * v[5159] + v[3449] * v[5160] + v[3450] * v[5161] - 2e0 * v[2329] * v[5458] * v[5459]) + (v[2166] * v[2250]
						+ v[2156] * v[2251] + v[2146] * v[2252] - v[8557 + i1663]) / 2e0 + v[3860] * v[935] + v[3862] * v[936] + v[3864] * v[937]
				+ v[3869] * v[944] + v[3871] * v[945] + v[3873] * v[950] - v[3855] * v[982] - v[3857] * v[983] + v[3856] * v[984]);
		v[5171] = -(v[2156] * v[3867]) + v[3961] - v[3866] * v[4724] - v[3319] * v[5297];
		v[3962] = v[3863] + v[3872];
		v[3963] = v[3861] + v[3865];
		v[3964] = v[3870] + v[3874];
		v[9758] = 0e0;
		v[9759] = 0e0;
		v[9760] = 0e0;
		v[9761] = 0e0;
		v[9762] = v[2206];
		v[9763] = v[2204];
		v[9764] = 0e0;
		v[9765] = 0e0;
		v[9766] = 0e0;
		v[9767] = 0e0;
		v[9768] = 0e0;
		v[9769] = 0e0;
		v[9770] = 0e0;
		v[9771] = 0e0;
		v[9772] = 0e0;
		v[9773] = 0e0;
		v[9774] = 0e0;
		v[9775] = 0e0;
		v[9776] = 0e0;
		v[9777] = 0e0;
		v[9778] = 0e0;
		v[9779] = 0e0;
		v[9780] = 0e0;
		v[9781] = 0e0;
		v[9710] = 0e0;
		v[9711] = 0e0;
		v[9712] = 0e0;
		v[9713] = v[2206];
		v[9714] = 0e0;
		v[9715] = v[2205];
		v[9716] = 0e0;
		v[9717] = 0e0;
		v[9718] = 0e0;
		v[9719] = 0e0;
		v[9720] = 0e0;
		v[9721] = 0e0;
		v[9722] = 0e0;
		v[9723] = 0e0;
		v[9724] = 0e0;
		v[9725] = 0e0;
		v[9726] = 0e0;
		v[9727] = 0e0;
		v[9728] = 0e0;
		v[9729] = 0e0;
		v[9730] = 0e0;
		v[9731] = 0e0;
		v[9732] = 0e0;
		v[9733] = 0e0;
		v[9686] = 0e0;
		v[9687] = 0e0;
		v[9688] = 0e0;
		v[9689] = v[2204];
		v[9690] = v[2205];
		v[9691] = 0e0;
		v[9692] = 0e0;
		v[9693] = 0e0;
		v[9694] = 0e0;
		v[9695] = 0e0;
		v[9696] = 0e0;
		v[9697] = 0e0;
		v[9698] = 0e0;
		v[9699] = 0e0;
		v[9700] = 0e0;
		v[9701] = 0e0;
		v[9702] = 0e0;
		v[9703] = 0e0;
		v[9704] = 0e0;
		v[9705] = 0e0;
		v[9706] = 0e0;
		v[9707] = 0e0;
		v[9708] = 0e0;
		v[9709] = 0e0;
		v[9614] = 0e0;
		v[9615] = 0e0;
		v[9616] = 0e0;
		v[9617] = 0e0;
		v[9618] = 0e0;
		v[9619] = 0e0;
		v[9620] = 0e0;
		v[9621] = 0e0;
		v[9622] = 0e0;
		v[9623] = 0e0;
		v[9624] = v[2203];
		v[9625] = v[2201];
		v[9626] = 0e0;
		v[9627] = 0e0;
		v[9628] = 0e0;
		v[9629] = 0e0;
		v[9630] = 0e0;
		v[9631] = 0e0;
		v[9632] = 0e0;
		v[9633] = 0e0;
		v[9634] = 0e0;
		v[9635] = 0e0;
		v[9636] = 0e0;
		v[9637] = 0e0;
		v[9566] = 0e0;
		v[9567] = 0e0;
		v[9568] = 0e0;
		v[9569] = 0e0;
		v[9570] = 0e0;
		v[9571] = 0e0;
		v[9572] = 0e0;
		v[9573] = 0e0;
		v[9574] = 0e0;
		v[9575] = v[2203];
		v[9576] = 0e0;
		v[9577] = v[2202];
		v[9578] = 0e0;
		v[9579] = 0e0;
		v[9580] = 0e0;
		v[9581] = 0e0;
		v[9582] = 0e0;
		v[9583] = 0e0;
		v[9584] = 0e0;
		v[9585] = 0e0;
		v[9586] = 0e0;
		v[9587] = 0e0;
		v[9588] = 0e0;
		v[9589] = 0e0;
		v[9542] = 0e0;
		v[9543] = 0e0;
		v[9544] = 0e0;
		v[9545] = 0e0;
		v[9546] = 0e0;
		v[9547] = 0e0;
		v[9548] = 0e0;
		v[9549] = 0e0;
		v[9550] = 0e0;
		v[9551] = v[2201];
		v[9552] = v[2202];
		v[9553] = 0e0;
		v[9554] = 0e0;
		v[9555] = 0e0;
		v[9556] = 0e0;
		v[9557] = 0e0;
		v[9558] = 0e0;
		v[9559] = 0e0;
		v[9560] = 0e0;
		v[9561] = 0e0;
		v[9562] = 0e0;
		v[9563] = 0e0;
		v[9564] = 0e0;
		v[9565] = 0e0;
		v[9470] = 0e0;
		v[9471] = 0e0;
		v[9472] = 0e0;
		v[9473] = 0e0;
		v[9474] = 0e0;
		v[9475] = 0e0;
		v[9476] = 0e0;
		v[9477] = 0e0;
		v[9478] = 0e0;
		v[9479] = 0e0;
		v[9480] = 0e0;
		v[9481] = 0e0;
		v[9482] = 0e0;
		v[9483] = 0e0;
		v[9484] = 0e0;
		v[9485] = 0e0;
		v[9486] = v[2200];
		v[9487] = v[2198];
		v[9488] = 0e0;
		v[9489] = 0e0;
		v[9490] = 0e0;
		v[9491] = 0e0;
		v[9492] = 0e0;
		v[9493] = 0e0;
		v[9422] = 0e0;
		v[9423] = 0e0;
		v[9424] = 0e0;
		v[9425] = 0e0;
		v[9426] = 0e0;
		v[9427] = 0e0;
		v[9428] = 0e0;
		v[9429] = 0e0;
		v[9430] = 0e0;
		v[9431] = 0e0;
		v[9432] = 0e0;
		v[9433] = 0e0;
		v[9434] = 0e0;
		v[9435] = 0e0;
		v[9436] = 0e0;
		v[9437] = v[2200];
		v[9438] = 0e0;
		v[9439] = v[2199];
		v[9440] = 0e0;
		v[9441] = 0e0;
		v[9442] = 0e0;
		v[9443] = 0e0;
		v[9444] = 0e0;
		v[9445] = 0e0;
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
		v[9410] = 0e0;
		v[9411] = 0e0;
		v[9412] = 0e0;
		v[9413] = v[2198];
		v[9414] = v[2199];
		v[9415] = 0e0;
		v[9416] = 0e0;
		v[9417] = 0e0;
		v[9418] = 0e0;
		v[9419] = 0e0;
		v[9420] = 0e0;
		v[9421] = 0e0;
		v[9326] = 0e0;
		v[9327] = 0e0;
		v[9328] = 0e0;
		v[9329] = 0e0;
		v[9330] = 0e0;
		v[9331] = 0e0;
		v[9332] = 0e0;
		v[9333] = 0e0;
		v[9334] = 0e0;
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
		v[9348] = v[2197];
		v[9349] = v[2195];
		v[9278] = 0e0;
		v[9279] = 0e0;
		v[9280] = 0e0;
		v[9281] = 0e0;
		v[9282] = 0e0;
		v[9283] = 0e0;
		v[9284] = 0e0;
		v[9285] = 0e0;
		v[9286] = 0e0;
		v[9287] = 0e0;
		v[9288] = 0e0;
		v[9289] = 0e0;
		v[9290] = 0e0;
		v[9291] = 0e0;
		v[9292] = 0e0;
		v[9293] = 0e0;
		v[9294] = 0e0;
		v[9295] = 0e0;
		v[9296] = 0e0;
		v[9297] = 0e0;
		v[9298] = 0e0;
		v[9299] = v[2197];
		v[9300] = 0e0;
		v[9301] = v[2196];
		v[9254] = 0e0;
		v[9255] = 0e0;
		v[9256] = 0e0;
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
		v[9273] = 0e0;
		v[9274] = 0e0;
		v[9275] = v[2195];
		v[9276] = v[2196];
		v[9277] = 0e0;
		v[9806] = v[3378] + v[19] * (-(v[1235] * v[2823]) - v[1211] * v[2824] + v[1283] * v[2825] + v[1259] * v[2829] + v[18] * (-
			(v[1566] * v[2836]) - v[1567] * v[2837] - v[1568] * v[2838]) + v[1409] * v[5027] - v[239] * v[5092]);
		v[9807] = v[3390] + v[19] * (-(v[1236] * v[2823]) - v[1212] * v[2824] + v[1284] * v[2825] + v[1260] * v[2829] + v[18] * (-
			(v[1570] * v[2836]) - v[1571] * v[2837] - v[1572] * v[2838]) + v[1409] * v[5026] - v[239] * v[5091]);
		v[9808] = v[3402] + v[19] * (-(v[1237] * v[2823]) - v[1213] * v[2824] + v[1285] * v[2825] + v[1261] * v[2829] + v[18] * (-
			(v[1574] * v[2836]) - v[1575] * v[2837] - v[1576] * v[2838]) + v[239] * v[5029]);
		v[9809] = v[3861] - v[3865] - v[3855] * v[4724] + 2e0 * (v[1407] * v[3205] + v[1403] * v[3211] + v[2288] * v[3975]
			+ v[3208] * v[4843]) + v[19] * (v[1747] * v[2527] + v[1746] * v[2531] + v[2329] * (v[3825] + v[1746] * v[3956]
				+ v[1747] * v[3958]) + v[3450] * v[3976] + v[3449] * v[4750] + v[3448] * v[4752] + v[2316] * v[5172] + v[2319] * v[5173]
				+ v[2322] * v[5174] + v[182] * v[5175] + v[187] * v[5176] + v[178] * v[5177]) + (v[1377] * v[3213] + v[1378] * v[3214]
					+ v[1379] * v[3315] - v[2181] * v[4972] + v[9757 + i1663]) / 2e0 + v[3962] * v[983] + v[3964] * v[984] + (-v[3859]
						+ v[5171]) * v[986];
		v[9810] = -v[3863] + v[3872] + v[3856] * v[4724] + 2e0 * (v[1407] * v[3203] + v[1402] * v[3214] + v[2287] * v[3977]
			+ v[3207] * v[4846]) + v[19] * (v[1748] * v[2527] + v[1746] * v[2533] + v[2209] * v[3448] + v[2210] * v[3450] + v[2329] *
				(v[3815] + v[1746] * v[3955] + v[1748] * v[3959]) + v[3449] * v[3978] + v[2323] * v[4948] + v[2317] * v[4950]
				+ v[2320] * v[5178] + v[179] * v[5179] + v[184] * v[5180] + v[2541] * v[5258]) + (v[1377] * v[3210] + v[1379] * v[3211]
					+ v[1378] * v[3316] + v[2182] * v[4972] + v[9709 + i1663]) / 2e0 + v[3964] * v[982] + v[3963] * v[983] + (-v[3859]
						- v[3876] + v[3961]) * v[985];
		v[9811] = v[3870] - v[3874] - v[3857] * v[4724] + 2e0 * (v[1403] * v[3210] + v[1402] * v[3213] + v[2286] * v[3979]
			+ v[3209] * v[4859]) + v[19] * (v[1748] * v[2531] + v[1747] * v[2533] + v[2211] * v[3449] + v[2212] * v[3450] + v[2329] *
				(v[3804] + v[1747] * v[3957] + v[1748] * v[3960]) + v[3448] * v[3980] + v[2321] * v[4949] + v[2318] * v[5181]
				+ v[2324] * v[5182] + v[190] * v[5183] + v[2545] * v[5256] + v[2549] * v[5257]) + (v[1378] * v[3203] + v[1379] * v[3205]
					+ v[1377] * v[3314] - v[2183] * v[4972] + v[9685 + i1663]) / 2e0 + (-v[3876] + v[5171]) * v[981] + v[3962] * v[982]
			+ v[3963] * v[984];
		v[9812] = v[3377] + v[19] * (-(v[1241] * v[2823]) - v[1217] * v[2824] + v[1289] * v[2825] + v[1265] * v[2829] + v[18] * (-
			(v[1590] * v[2836]) - v[1591] * v[2837] - v[1592] * v[2838]) + v[1410] * v[5027] - v[240] * v[5092]);
		v[9813] = v[3389] + v[19] * (-(v[1242] * v[2823]) - v[1218] * v[2824] + v[1290] * v[2825] + v[1266] * v[2829] + v[18] * (-
			(v[1594] * v[2836]) - v[1595] * v[2837] - v[1596] * v[2838]) + v[1410] * v[5026] - v[240] * v[5091]);
		v[9814] = v[3401] + v[19] * (-(v[1243] * v[2823]) - v[1219] * v[2824] + v[1291] * v[2825] + v[1267] * v[2829] + v[18] * (-
			(v[1598] * v[2836]) - v[1599] * v[2837] - v[1600] * v[2838]) + v[240] * v[5029]);
		v[9815] = v[3839] - v[3843] - v[3833] * v[4728] + 2e0 * (v[1401] * v[3216] + v[1397] * v[3221] + v[2285] * v[3981]
			+ v[3218] * v[4837]) + v[19] * (v[1744] * v[2555] + v[1743] * v[2559] + v[2351] * (v[3795] + v[1743] * v[3934]
				+ v[1744] * v[3936]) + v[3447] * v[3982] + v[3446] * v[4755] + v[3445] * v[4757] + v[2338] * v[5185] + v[2341] * v[5186]
				+ v[2344] * v[5187] + v[201] * v[5188] + v[206] * v[5189] + v[197] * v[5190]) + (v[1371] * v[3222] + v[1372] * v[3223]
					+ v[1373] * v[3327] - v[2167] * v[4976] + v[9613 + i1663]) / 2e0 + v[3940] * v[977] + v[3942] * v[978] + (-v[3837]
						+ v[5184]) * v[980];
		v[9816] = -v[3841] + v[3850] + v[3834] * v[4728] + 2e0 * (v[1401] * v[3215] + v[1396] * v[3223] + v[2284] * v[3983]
			+ v[3217] * v[4838]) + v[19] * (v[1745] * v[2555] + v[1743] * v[2561] + v[2215] * v[3445] + v[2216] * v[3447] + v[2351] *
				(v[3785] + v[1743] * v[3933] + v[1745] * v[3937]) + v[3446] * v[3984] + v[2345] * v[4945] + v[2339] * v[4947]
				+ v[2342] * v[5191] + v[198] * v[5192] + v[203] * v[5193] + v[2569] * v[5255]) + (v[1371] * v[3220] + v[1373] * v[3221]
					+ v[1372] * v[3328] + v[2168] * v[4976] + v[9565 + i1663]) / 2e0 + v[3942] * v[976] + v[3941] * v[977] + (-v[3837]
						- v[3854] + v[3939]) * v[979];
		v[9817] = v[3848] - v[3852] - v[3835] * v[4728] + 2e0 * (v[1397] * v[3220] + v[1396] * v[3222] + v[2283] * v[3985]
			+ v[3219] * v[4839]) + v[19] * (v[1745] * v[2559] + v[1744] * v[2561] + v[2217] * v[3446] + v[2218] * v[3447] + v[2351] *
				(v[3774] + v[1744] * v[3935] + v[1745] * v[3938]) + v[3445] * v[3986] + v[2343] * v[4946] + v[2340] * v[5194]
				+ v[2346] * v[5195] + v[209] * v[5196] + v[2573] * v[5253] + v[2577] * v[5254]) + (v[1372] * v[3215] + v[1373] * v[3216]
					+ v[1371] * v[3326] - v[2169] * v[4976] + v[9541 + i1663]) / 2e0 + (-v[3854] + v[5184]) * v[975] + v[3940] * v[976]
			+ v[3941] * v[978];
		v[9818] = v[3372] + v[19] * (-(v[1247] * v[2823]) - v[1223] * v[2824] + v[1295] * v[2825] + v[1271] * v[2829] + v[18] * (-
			(v[1614] * v[2836]) - v[1615] * v[2837] - v[1616] * v[2838]) + v[1411] * v[5027] + v[326] * v[5092]);
		v[9819] = v[3384] + v[19] * (-(v[1248] * v[2823]) - v[1224] * v[2824] + v[1296] * v[2825] + v[1272] * v[2829] + v[18] * (-
			(v[1618] * v[2836]) - v[1619] * v[2837] - v[1620] * v[2838]) + v[1411] * v[5026] + v[326] * v[5091]);
		v[9820] = v[3396] + v[19] * (-(v[1249] * v[2823]) - v[1225] * v[2824] + v[1297] * v[2825] + v[1273] * v[2829] + v[18] * (-
			(v[1622] * v[2836]) - v[1623] * v[2837] - v[1624] * v[2838]) - v[326] * v[5029]);
		v[9821] = v[3727] - v[3731] - v[3721] * v[4737] + 2e0 * (v[1395] * v[3225] + v[1391] * v[3230] + v[2282] * v[3987]
			+ v[3227] * v[4831]) + v[19] * (v[1741] * v[2583] + v[1740] * v[2587] + v[2433] * (v[3691] + v[1740] * v[3912]
				+ v[1741] * v[3914]) + v[3444] * v[3988] + v[3443] * v[4760] + v[3442] * v[4762] + v[2420] * v[5198] + v[2423] * v[5199]
				+ v[2426] * v[5200] + v[269] * v[5201] + v[274] * v[5202] + v[265] * v[5203]) + (v[1365] * v[3231] + v[1366] * v[3232]
					+ v[1367] * v[3339] - v[2084] * v[4986] + v[9469 + i1663]) / 2e0 + v[3918] * v[971] + v[3920] * v[972] + (-v[3725]
						+ v[5197]) * v[974];
		v[9822] = -v[3729] + v[3738] + v[3722] * v[4737] + 2e0 * (v[1395] * v[3224] + v[1390] * v[3232] + v[2281] * v[3989]
			+ v[3226] * v[4832]) + v[19] * (v[1742] * v[2583] + v[1740] * v[2589] + v[2221] * v[3442] + v[2222] * v[3444] + v[2433] *
				(v[3681] + v[1740] * v[3911] + v[1742] * v[3915]) + v[3443] * v[3990] + v[2427] * v[4916] + v[2421] * v[4918]
				+ v[2424] * v[5204] + v[266] * v[5205] + v[271] * v[5206] + v[2597] * v[5252]) + (v[1365] * v[3229] + v[1367] * v[3230]
					+ v[1366] * v[3340] + v[2085] * v[4986] + v[9421 + i1663]) / 2e0 + v[3920] * v[970] + v[3919] * v[971] + (-v[3725]
						- v[3742] + v[3917]) * v[973];
		v[9823] = v[3736] - v[3740] - v[3723] * v[4737] + 2e0 * (v[1391] * v[3229] + v[1390] * v[3231] + v[2280] * v[3991]
			+ v[3228] * v[4833]) + v[19] * (v[1742] * v[2587] + v[1741] * v[2589] + v[2223] * v[3443] + v[2224] * v[3444] + v[2433] *
				(v[3670] + v[1741] * v[3913] + v[1742] * v[3916]) + v[3442] * v[3992] + v[2425] * v[4917] + v[2422] * v[5207]
				+ v[2428] * v[5208] + v[277] * v[5209] + v[2601] * v[5250] + v[2605] * v[5251]) + (v[1366] * v[3224] + v[1367] * v[3225]
					+ v[1365] * v[3338] - v[2086] * v[4986] + v[9397 + i1663]) / 2e0 + (-v[3742] + v[5197]) * v[969] + v[3918] * v[970]
			+ v[3919] * v[972];
		v[9824] = v[3371] + v[19] * (-(v[1253] * v[2823]) - v[1229] * v[2824] + v[1301] * v[2825] + v[1277] * v[2829] + v[18] * (-
			(v[1638] * v[2836]) - v[1639] * v[2837] - v[1640] * v[2838]) + v[1412] * v[5027] + v[327] * v[5092]);
		v[9825] = v[3383] + v[19] * (-(v[1254] * v[2823]) - v[1230] * v[2824] + v[1302] * v[2825] + v[1278] * v[2829] + v[18] * (-
			(v[1642] * v[2836]) - v[1643] * v[2837] - v[1644] * v[2838]) + v[1412] * v[5026] + v[327] * v[5091]);
		v[9826] = v[3395] + v[19] * (-(v[1255] * v[2823]) - v[1231] * v[2824] + v[1303] * v[2825] + v[1279] * v[2829] + v[18] * (-
			(v[1646] * v[2836]) - v[1647] * v[2837] - v[1648] * v[2838]) - v[327] * v[5029]);
		v[9827] = v[3705] - v[3709] - v[3699] * v[4741] + 2e0 * (v[1389] * v[3234] + v[1385] * v[3239] + v[2279] * v[3993]
			+ v[3236] * v[4825]) + v[19] * (v[1738] * v[2611] + v[1737] * v[2615] + v[2455] * (v[3661] + v[1737] * v[3890]
				+ v[1738] * v[3892]) + v[3441] * v[3994] + v[3440] * v[4765] + v[3439] * v[4767] + v[2442] * v[5211] + v[2445] * v[5212]
				+ v[2448] * v[5213] + v[288] * v[5214] + v[293] * v[5215] + v[284] * v[5216]) + (v[1359] * v[3240] + v[1360] * v[3241]
					+ v[1361] * v[3351] - v[2070] * v[4990] + v[9325 + i1663]) / 2e0 + v[3896] * v[965] + v[3898] * v[966] + (-v[3703]
						+ v[5210]) * v[968];
		v[9828] = -v[3707] + v[3716] + v[3700] * v[4741] + 2e0 * (v[1389] * v[3233] + v[1384] * v[3241] + v[2278] * v[3995]
			+ v[3235] * v[4826]) + v[19] * (v[1739] * v[2611] + v[1737] * v[2617] + v[2227] * v[3439] + v[2228] * v[3441] + v[2455] *
				(v[3651] + v[1737] * v[3889] + v[1739] * v[3893]) + v[3440] * v[3996] + v[2449] * v[4913] + v[2443] * v[4915]
				+ v[2446] * v[5217] + v[285] * v[5218] + v[290] * v[5219] + v[2625] * v[5249]) + (v[1359] * v[3238] + v[1361] * v[3239]
					+ v[1360] * v[3352] + v[2071] * v[4990] + v[9277 + i1663]) / 2e0 + v[3898] * v[964] + v[3897] * v[965] + (-v[3703]
						- v[3720] + v[3895]) * v[967];
		v[9829] = v[3714] - v[3718] - v[3701] * v[4741] + 2e0 * (v[1385] * v[3238] + v[1384] * v[3240] + v[2277] * v[3997]
			+ v[3237] * v[4827]) + v[19] * (v[1739] * v[2615] + v[1738] * v[2617] + v[2229] * v[3440] + v[2230] * v[3441] + v[2455] *
				(v[3640] + v[1738] * v[3891] + v[1739] * v[3894]) + v[3439] * v[3998] + v[2447] * v[4914] + v[2444] * v[5220]
				+ v[2450] * v[5221] + v[296] * v[5222] + v[2629] * v[5247] + v[2633] * v[5248]) + (v[1360] * v[3233] + v[1361] * v[3234]
					+ v[1359] * v[3350] - v[2072] * v[4990] + v[9253 + i1663]) / 2e0 + (-v[3720] + v[5210]) * v[963] + v[3896] * v[964]
			+ v[3897] * v[966];
		v[3965] = -(v[3622] * v[3623]);
		v[3969] = -(v[1462] * (v[335] * (v[139] * v[3619] + v[5166]) + v[1465] * (v[3965] * v[4699] + v[3968] * v[5167])
			+ v[3965] * v[5226] + v[5165] * v[535])) + v[1565] * (v[335] * (v[140] * v[3618] - v[5165]) + v[1463] *
				(v[3965] * v[4697] + v[3966] * v[5167]) + v[3965] * v[5223] + v[5166] * v[535]);
		v[3970] = -(v[3756] * v[3757]);
		v[3974] = -(v[1457] * (v[507] * v[5168] + v[1460] * (v[3970] * v[4705] + v[3973] * v[5170]) + v[3970] * v[5230]
			+ v[248] * (v[5169] + v[3753] * v[76]))) + v[1564] * (v[507] * v[5169] + v[1458] * (v[3970] * v[4703] + v[3971] * v[5170]
				) + v[3970] * v[5227] + v[248] * (-v[5168] + v[3752] * v[77]));
		Rc[i1663 - 1] += v[6126 + i1663];
		for (i2234 = 1; i2234 <= 24; i2234++) {
			Kc[i1663 - 1][i2234 - 1] += v[9805 + i2234] + v[3974] * v[9829 + i2234] + v[3969] * v[9853 + i2234] + v[3766] * v[9877
				+ i2234] + v[3632] * v[9901 + i2234];
		};/* end for */
	};/* end for */
	v[4004] = v[1347] * v[1485];
	v[4005] = v[1348] * v[1485];
	v[4006] = v[1349] * v[1485];
	v[4007] = v[1350] * v[1485];
	v[4008] = v[1351] * v[1485];
	v[4009] = v[1352] * v[1485];
	v[4010] = v[1353] * v[1485];
	v[4011] = v[1354] * v[1485];
	v[4012] = v[1355] * v[1485];
	v[4013] = v[1356] * v[1485];
	v[4014] = v[1357] * v[1485];
	v[4015] = v[1358] * v[1485];
	v[4016] = v[1347] * v[1484];
	v[4017] = v[1348] * v[1484];
	v[4018] = v[1349] * v[1484];
	v[4019] = v[1350] * v[1484];
	v[4020] = v[1351] * v[1484];
	v[4021] = v[1352] * v[1484];
	v[4022] = v[1353] * v[1484];
	v[4023] = v[1354] * v[1484];
	v[4024] = v[1355] * v[1484];
	v[4025] = v[1356] * v[1484];
	v[4026] = v[1357] * v[1484];
	v[4027] = v[1358] * v[1484];
	v[4031] = -(v[1483] * v[4028]) - v[1484] * v[4029] - v[1485] * v[4030];
	v[4035] = -(v[1483] * v[4032]) - v[1484] * v[4033] - v[1485] * v[4034];
	v[5464] = v[140] * v[1565] * v[4031] - v[139] * v[1462] * v[4035];
	v[4036] = v[1347] * v[1483];
	v[4037] = v[1348] * v[1483];
	v[4038] = v[1349] * v[1483];
	v[4673] = v[4004] * v[4784] + v[4017] * v[4785] + v[4038] * v[4786] + v[4005] * v[722] + v[4006] * v[723]
		+ v[4016] * v[724] + v[4018] * v[731] + v[4036] * v[732] + v[4037] * v[737];
	v[4039] = v[1350] * v[1483];
	v[4040] = v[1351] * v[1483];
	v[4041] = v[1352] * v[1483];
	v[4679] = v[4007] * v[4787] + v[4020] * v[4788] + v[4041] * v[4789] + v[4008] * v[753] + v[4009] * v[754]
		+ v[4019] * v[755] + v[4021] * v[762] + v[4039] * v[763] + v[4040] * v[768];
	v[4045] = v[1483] * v[4042] + v[1484] * v[4043] + v[1485] * v[4044];
	v[4049] = v[1483] * v[4046] + v[1484] * v[4047] + v[1485] * v[4048];
	v[4050] = v[1353] * v[1483];
	v[4051] = v[1354] * v[1483];
	v[4052] = v[1355] * v[1483];
	v[4685] = v[4010] * v[4796] + v[4023] * v[4797] + v[4052] * v[4798] + v[4011] * v[904] + v[4012] * v[905]
		+ v[4022] * v[906] + v[4024] * v[913] + v[4050] * v[914] + v[4051] * v[919];
	v[4053] = v[1356] * v[1483];
	v[4054] = v[1357] * v[1483];
	v[4055] = v[1358] * v[1483];
	v[4691] = v[4013] * v[4799] + v[4026] * v[4800] + v[4055] * v[4801] + v[4014] * v[935] + v[4015] * v[936]
		+ v[4025] * v[937] + v[4027] * v[944] + v[4053] * v[945] + v[4054] * v[950];
	v[4056] = v[4031] * v[4748] + v[4035] * v[5087];
	v[5225] = v[3623] * v[4056];
	v[4700] = v[4035] * v[4747] + v[5223] * v[5225];
	v[4696] = v[4031] * v[5224] - v[5225] * v[5226];
	v[4057] = (v[1483] * v[520] - v[1483] * v[521] + v[1484] * v[523] - v[1484] * v[524] + v[1485] * v[526] - v[1485] * v[527])
		/ 2e0;
	v[4059] = v[281] * v[4005];
	v[4060] = v[281] * v[4006];
	v[4061] = v[281] * v[4016];
	v[4062] = v[4017] * v[4741];
	v[4063] = v[281] * v[4018];
	v[4064] = v[281] * v[4036];
	v[4065] = v[281] * v[4037];
	v[4066] = v[4038] * v[4741];
	v[4719] = -v[4062] - v[4066] + v[4673] * v[749];
	v[4718] = v[4062] + v[4719] - v[4004] * v[4741];
	v[4717] = -v[4062] + v[4066] + v[4718];
	v[4069] = v[262] * v[4008];
	v[4070] = v[262] * v[4009];
	v[4071] = v[262] * v[4019];
	v[4072] = v[4020] * v[4737];
	v[4073] = v[262] * v[4021];
	v[4074] = v[262] * v[4039];
	v[4075] = v[262] * v[4040];
	v[4076] = v[4041] * v[4737];
	v[4716] = -v[4072] - v[4076] + v[4679] * v[780];
	v[4715] = v[4072] + v[4716] - v[4007] * v[4737];
	v[4714] = -v[4072] + v[4076] + v[4715];
	v[4078] = v[4045] * v[4735] + v[4049] * v[5105];
	v[5229] = v[3757] * v[4078];
	v[4706] = v[4049] * v[4734] + v[5227] * v[5229];
	v[4702] = v[4045] * v[5228] - v[5229] * v[5230];
	v[4079] = (-(v[1483] * v[492]) + v[1483] * v[493] - v[1484] * v[495] + v[1484] * v[496] - v[1485] * v[498]
		+ v[1485] * v[499]) / 2e0;
	v[4081] = v[194] * v[4011];
	v[4082] = v[194] * v[4012];
	v[4083] = v[194] * v[4022];
	v[4084] = v[4023] * v[4728];
	v[4085] = v[194] * v[4024];
	v[4086] = v[194] * v[4050];
	v[4087] = v[194] * v[4051];
	v[4088] = v[4052] * v[4728];
	v[4713] = -v[4084] - v[4088] + v[4685] * v[931];
	v[4712] = v[4084] + v[4713] - v[4010] * v[4728];
	v[4711] = -v[4084] + v[4088] + v[4712];
	v[4091] = v[175] * v[4014];
	v[4092] = v[175] * v[4015];
	v[4093] = v[175] * v[4025];
	v[4094] = v[4026] * v[4724];
	v[4095] = v[175] * v[4027];
	v[4096] = v[175] * v[4053];
	v[4097] = v[175] * v[4054];
	v[4098] = v[4055] * v[4724];
	v[4710] = -v[4094] - v[4098] + v[4691] * v[962];
	v[4709] = v[4094] + v[4710] - v[4013] * v[4724];
	v[4708] = -v[4094] + v[4098] + v[4709];
	v[4100] = v[4060] + v[4064];
	v[4101] = v[4059] + v[4061];
	v[4102] = v[4063] + v[4065];
	v[4103] = v[4070] + v[4074];
	v[4104] = v[4069] + v[4071];
	v[4105] = v[4073] + v[4075];
	v[4106] = v[4082] + v[4086];
	v[4107] = v[4081] + v[4083];
	v[4108] = v[4085] + v[4087];
	v[4109] = v[4092] + v[4096];
	v[4110] = v[4091] + v[4093];
	v[4111] = v[4095] + v[4097];
	v[9933] = v[5280];
	v[9934] = v[5282];
	v[9935] = v[5284];
	v[9936] = v[4091] - v[4093] + v[4109] * v[983] + v[4111] * v[984] + v[4708] * v[986];
	v[9937] = -v[4092] + v[4096] + v[4111] * v[982] + v[4110] * v[983] + v[4709] * v[985];
	v[9938] = v[4095] - v[4097] + v[4710] * v[981] + v[4109] * v[982] + v[4110] * v[984];
	v[9939] = v[5281];
	v[9940] = v[5283];
	v[9941] = v[5285];
	v[9942] = v[4081] - v[4083] + v[4106] * v[977] + v[4108] * v[978] + v[4711] * v[980];
	v[9943] = -v[4082] + v[4086] + v[4108] * v[976] + v[4107] * v[977] + v[4712] * v[979];
	v[9944] = v[4085] - v[4087] + v[4713] * v[975] + v[4106] * v[976] + v[4107] * v[978];
	v[9945] = v[5259];
	v[9946] = v[5261];
	v[9947] = v[5263];
	v[9948] = v[4069] - v[4071] + v[4103] * v[971] + v[4105] * v[972] + v[4714] * v[974];
	v[9949] = -v[4070] + v[4074] + v[4105] * v[970] + v[4104] * v[971] + v[4715] * v[973];
	v[9950] = v[4073] - v[4075] + v[4716] * v[969] + v[4103] * v[970] + v[4104] * v[972];
	v[9951] = v[5260];
	v[9952] = v[5262];
	v[9953] = v[5264];
	v[9954] = v[4059] - v[4061] + v[4100] * v[965] + v[4102] * v[966] + v[4717] * v[968];
	v[9955] = -v[4060] + v[4064] + v[4102] * v[964] + v[4101] * v[965] + v[4718] * v[967];
	v[9956] = v[4063] - v[4065] + v[4719] * v[963] + v[4100] * v[964] + v[4101] * v[966];
	v[4113] = v[1565] * v[4696] - v[1462] * v[4700];
	v[4115] = v[1564] * v[4702] - v[1457] * v[4706];
	for (i4002 = 1; i4002 <= 24; i4002++) {
		v[4265] = v[9901 + i4002];
		v[5231] = v[4265] / 2e0;
		v[4433] = v[1483] * v[5231];
		v[4427] = v[1484] * v[5231];
		v[4421] = v[1485] * v[5231];
		v[4264] = v[9877 + i4002];
		v[5232] = -0.5e0 * v[4264];
		v[4436] = v[1483] * v[5232];
		v[4430] = v[1484] * v[5232];
		v[4424] = v[1485] * v[5232];
		v[4218] = (i4002 == 23 ? 1 : 0);
		v[4215] = (i4002 == 22 ? 1 : 0);
		v[4212] = (i4002 == 24 ? 1 : 0);
		v[4209] = (i4002 == 17 ? 1 : 0);
		v[4206] = (i4002 == 16 ? 1 : 0);
		v[4203] = (i4002 == 18 ? 1 : 0);
		v[4200] = (i4002 == 11 ? 1 : 0);
		v[4197] = (i4002 == 10 ? 1 : 0);
		v[4194] = (i4002 == 12 ? 1 : 0);
		v[4191] = (i4002 == 5 ? 1 : 0);
		v[4188] = (i4002 == 4 ? 1 : 0);
		v[4185] = (i4002 == 6 ? 1 : 0);
		v[4182] = v[9853 + i4002];
		v[4164] = v[9829 + i4002];
		v[4125] = v[6157 + i4002];
		v[4127] = v[6181 + i4002];
		v[4129] = v[6205 + i4002];
		v[4131] = v[6229 + i4002];
		v[4133] = v[6253 + i4002];
		v[4135] = v[6277 + i4002];
		v[4137] = v[6301 + i4002];
		v[4139] = v[6325 + i4002];
		v[4141] = v[6349 + i4002];
		v[4143] = v[6373 + i4002];
		v[4145] = v[6397 + i4002];
		v[4147] = v[6421 + i4002];
		v[4149] = v[6445 + i4002];
		v[4666] = (v[4149] * v[962]) / 2e0;
		v[5233] = 2e0 * v[4666];
		v[4151] = v[6469 + i4002];
		v[4153] = v[6565 + i4002];
		v[4155] = v[6661 + i4002];
		v[4157] = v[6757 + i4002];
		v[4653] = (v[4157] * v[931]) / 2e0;
		v[5234] = 2e0 * v[4653];
		v[4159] = v[6781 + i4002];
		v[4161] = v[6877 + i4002];
		v[4163] = v[6973 + i4002];
		v[4165] = v[4164] * v[507];
		v[4167] = v[7213 + i4002];
		v[4585] = (v[4167] * v[780]) / 2e0;
		v[5235] = 2e0 * v[4585];
		v[4169] = v[7237 + i4002];
		v[4171] = v[7333 + i4002];
		v[4173] = v[7429 + i4002];
		v[4175] = v[7525 + i4002];
		v[4572] = (v[4175] * v[749]) / 2e0;
		v[5236] = 2e0 * v[4572];
		v[4177] = v[7549 + i4002];
		v[4179] = v[7645 + i4002];
		v[4181] = v[7741 + i4002];
		v[4183] = v[4182] * v[535];
		v[4184] = v[4125] - v[4185];
		v[4186] = v[4125] + v[4185];
		v[4187] = v[4127] - v[4188];
		v[4189] = v[4127] + v[4188];
		v[4190] = v[4129] + v[4191];
		v[4192] = v[4129] - v[4191];
		v[4193] = v[4131] - v[4194];
		v[4195] = v[4131] + v[4194];
		v[4196] = v[4133] - v[4197];
		v[4198] = v[4133] + v[4197];
		v[4199] = v[4135] + v[4200];
		v[4201] = v[4135] - v[4200];
		v[4202] = v[4137] - v[4203];
		v[4204] = v[4137] + v[4203];
		v[4205] = v[4139] - v[4206];
		v[4207] = v[4139] + v[4206];
		v[4208] = v[4141] + v[4209];
		v[4210] = v[4141] - v[4209];
		v[4211] = v[4143] - v[4212];
		v[4213] = v[4143] + v[4212];
		v[4214] = v[4145] - v[4215];
		v[4216] = v[4145] + v[4215];
		v[4217] = v[4147] + v[4218];
		v[4219] = v[4147] - v[4218];
		v[4220] = v[4151] * v[4724] + v[4801] * v[5233];
		v[4221] = v[175] * v[4184] + v[5233] * v[950];
		v[4222] = v[175] * v[4190] + v[5233] * v[945];
		v[4223] = v[175] * v[4186] + v[5233] * v[944];
		v[4224] = v[4153] * v[4724] + v[4800] * v[5233];
		v[4225] = v[175] * v[4187] + v[5233] * v[937];
		v[4226] = v[175] * v[4192] + v[5233] * v[936];
		v[4227] = v[175] * v[4189] + v[5233] * v[935];
		v[4228] = v[4155] * v[4724] + v[4799] * v[5233];
		v[4229] = v[4159] * v[4728] + v[4798] * v[5234];
		v[4230] = v[194] * v[4193] + v[5234] * v[919];
		v[4231] = v[194] * v[4199] + v[5234] * v[914];
		v[4232] = v[194] * v[4195] + v[5234] * v[913];
		v[4233] = v[4161] * v[4728] + v[4797] * v[5234];
		v[4234] = v[194] * v[4196] + v[5234] * v[906];
		v[4235] = v[194] * v[4201] + v[5234] * v[905];
		v[4236] = v[194] * v[4198] + v[5234] * v[904];
		v[4237] = v[4163] * v[4728] + v[4796] * v[5234];
		v[4239] = v[4164] * v[4238] + v[4165] * v[5105];
		v[4241] = v[4164] * v[4240] + v[4165] * v[4735];
		v[4242] = v[4169] * v[4737] + v[4789] * v[5235];
		v[4243] = v[262] * v[4202] + v[5235] * v[768];
		v[4244] = v[262] * v[4208] + v[5235] * v[763];
		v[4245] = v[262] * v[4204] + v[5235] * v[762];
		v[4246] = v[4171] * v[4737] + v[4788] * v[5235];
		v[4247] = v[262] * v[4205] + v[5235] * v[755];
		v[4248] = v[262] * v[4210] + v[5235] * v[754];
		v[4249] = v[262] * v[4207] + v[5235] * v[753];
		v[4250] = v[4173] * v[4737] + v[4787] * v[5235];
		v[4251] = v[4177] * v[4741] + v[4786] * v[5236];
		v[4252] = v[281] * v[4211] + v[5236] * v[737];
		v[4253] = v[281] * v[4217] + v[5236] * v[732];
		v[4254] = v[281] * v[4213] + v[5236] * v[731];
		v[4255] = v[4179] * v[4741] + v[4785] * v[5236];
		v[4256] = v[281] * v[4214] + v[5236] * v[724];
		v[4257] = v[281] * v[4219] + v[5236] * v[723];
		v[4258] = v[281] * v[4216] + v[5236] * v[722];
		v[4259] = v[4181] * v[4741] + v[4784] * v[5236];
		v[4261] = v[4182] * v[4260] + v[4183] * v[5087];
		v[4263] = v[4182] * v[4262] + v[4183] * v[4748];
		v[4266] = v[10035 + i4002] + v[1358] * v[4220] + v[1357] * v[4221] + v[1356] * v[4222] + v[1355] * v[4229]
			+ v[1354] * v[4230] + v[1353] * v[4231] + v[4046] * v[4239] + v[4042] * v[4241] + v[1352] * v[4242] + v[1351] * v[4243]
			+ v[1350] * v[4244] + v[1349] * v[4251] + v[1348] * v[4252] + v[1347] * v[4253] - v[4032] * v[4261] - v[4028] * v[4263]
			+ v[4264] * v[494] - v[4265] * v[522];
		v[5237] = v[17] * v[4266];
		v[4397] = v[480] * v[5237];
		v[4282] = v[1383] * v[5237];
		v[4267] = v[10011 + i4002] + v[1358] * v[4223] + v[1357] * v[4224] + v[1356] * v[4225] + v[1355] * v[4232]
			+ v[1354] * v[4233] + v[1353] * v[4234] + v[4047] * v[4239] + v[4043] * v[4241] + v[1352] * v[4245] + v[1351] * v[4246]
			+ v[1350] * v[4247] + v[1349] * v[4254] + v[1348] * v[4255] + v[1347] * v[4256] - v[4033] * v[4261] - v[4029] * v[4263]
			+ v[4264] * v[497] - v[4265] * v[525];
		v[5243] = v[4267] * v[4275] + v[4282];
		v[5239] = v[4266] * v[480] + v[4267] * v[481];
		v[5238] = v[17] * v[4267];
		v[4392] = v[481] * v[5238];
		v[4279] = v[1413] * v[5238];
		v[5241] = v[4266] * v[4275] + v[4279];
		v[4268] = v[1358] * v[4226] + v[1357] * v[4227] + v[1356] * v[4228] + v[1355] * v[4235] + v[1354] * v[4236]
			+ v[1353] * v[4237] + v[4048] * v[4239] + v[4044] * v[4241] + v[1352] * v[4248] + v[1351] * v[4249] + v[1350] * v[4250]
			+ v[1349] * v[4257] + v[1348] * v[4258] + v[1347] * v[4259] - v[4034] * v[4261] - v[4030] * v[4263] + v[4264] * v[500]
			- v[4265] * v[528] + v[9987 + i4002];
		v[5267] = v[17] * v[4268];
		v[5242] = v[4268] * v[480];
		v[5240] = v[4268] * v[481];
		v[4720] = v[4268] * v[4773];
		v[5322] = v[4279] + v[4720] * v[481];
		v[5321] = v[4282] + v[4720] * v[480];
		v[4272] = v[1455] * v[5267];
		v[4270] = v[376] * v[4272] + v[4269] * v[5239];
		v[4273] = v[366] * v[4272] + v[4271] * v[5239];
		v[4274] = v[4269] * v[5240] + v[376] * v[5241];
		v[4276] = v[4271] * v[5242] + v[366] * v[5243];
		v[4278] = v[386] * v[4272] + v[4277] * v[5239];
		v[4280] = v[4277] * v[5240] + v[386] * v[5241];
		v[4281] = v[4271] * v[5240] + v[366] * v[5241];
		v[4283] = v[4277] * v[5242] + v[386] * v[5243];
		v[4284] = v[4269] * v[5242] + v[376] * v[5243];
		v[4286] = v[402] * v[4272] + v[4285] * v[5239];
		v[4288] = v[392] * v[4272] + v[4287] * v[5239];
		v[4289] = v[4285] * v[5240] + v[402] * v[5241];
		v[4290] = v[4287] * v[5242] + v[392] * v[5243];
		v[4292] = v[412] * v[4272] + v[4291] * v[5239];
		v[4293] = v[4291] * v[5240] + v[412] * v[5241];
		v[4294] = v[4287] * v[5240] + v[392] * v[5241];
		v[4295] = v[4291] * v[5242] + v[412] * v[5243];
		v[4296] = v[4285] * v[5242] + v[402] * v[5243];
		v[4298] = v[4272] * v[428] + v[4297] * v[5239];
		v[4300] = v[418] * v[4272] + v[4299] * v[5239];
		v[4301] = v[4297] * v[5240] + v[428] * v[5241];
		v[4302] = v[4299] * v[5242] + v[418] * v[5243];
		v[4304] = v[4272] * v[438] + v[4303] * v[5239];
		v[4305] = v[4303] * v[5240] + v[438] * v[5241];
		v[4306] = v[4299] * v[5240] + v[418] * v[5241];
		v[4307] = v[4303] * v[5242] + v[438] * v[5243];
		v[4308] = v[4297] * v[5242] + v[428] * v[5243];
		v[4310] = v[4272] * v[454] + v[4309] * v[5239];
		v[4312] = v[4272] * v[444] + v[4311] * v[5239];
		v[4313] = v[4309] * v[5240] + v[454] * v[5241];
		v[4314] = v[4311] * v[5242] + v[444] * v[5243];
		v[4316] = v[4272] * v[464] + v[4315] * v[5239];
		v[4317] = v[4315] * v[5240] + v[464] * v[5241];
		v[4318] = v[4311] * v[5240] + v[444] * v[5241];
		v[4319] = v[4315] * v[5242] + v[464] * v[5243];
		v[4320] = v[4309] * v[5242] + v[454] * v[5243];
		v[4321] = v[4392] + v[4397];
		v[5266] = v[4321] * v[482];
		v[5323] = v[4272] + v[5266];
		v[4322] = v[17] * (v[350] * v[4266] + v[349] * v[4267]);
		v[4323] = v[17] * (v[353] * v[4266] + v[352] * v[4267]);
		v[4324] = v[17] * (v[356] * v[4266] + v[355] * v[4267]);
		v[4325] = v[17] * (v[359] * v[4266] + v[358] * v[4267]);
		v[5244] = v[239] * v[4322] + v[240] * v[4323] - v[326] * v[4324] - v[327] * v[4325];
		v[4326] = v[17] * (v[3277] * v[4266] + v[3167] * v[4267] + v[3065] * v[4268]);
		v[4544] = v[2230] * v[4326];
		v[4530] = v[3998] * v[4326];
		v[4327] = v[17] * (v[3279] * v[4266] + v[3169] * v[4267] + v[3067] * v[4268]);
		v[4537] = v[3996] * v[4327];
		v[5269] = v[2229] * v[4326] + v[4537];
		v[5268] = v[2227] * v[4327] + v[4530];
		v[4328] = v[17] * (v[3281] * v[4266] + v[3171] * v[4267] + v[3069] * v[4268]);
		v[4542] = v[3994] * v[4328];
		v[5270] = v[2228] * v[4327] + v[4542];
		v[4539] = v[4328] * v[4765];
		v[4533] = v[4328] * v[4767];
		v[4329] = v[17] * (v[3283] * v[4266] + v[3173] * v[4267] + v[3071] * v[4268]);
		v[4562] = v[2224] * v[4329];
		v[4548] = v[3992] * v[4329];
		v[4330] = v[17] * (v[3285] * v[4266] + v[3175] * v[4267] + v[3073] * v[4268]);
		v[4555] = v[3990] * v[4330];
		v[5272] = v[2223] * v[4329] + v[4555];
		v[5271] = v[2221] * v[4330] + v[4548];
		v[4331] = v[17] * (v[3287] * v[4266] + v[3177] * v[4267] + v[3075] * v[4268]);
		v[4560] = v[3988] * v[4331];
		v[5273] = v[2222] * v[4330] + v[4560];
		v[4557] = v[4331] * v[4760];
		v[4551] = v[4331] * v[4762];
		v[4332] = v[17] * (v[3289] * v[4266] + v[3179] * v[4267] + v[3077] * v[4268]);
		v[4625] = v[2218] * v[4332];
		v[4611] = v[3986] * v[4332];
		v[4333] = v[17] * (v[3291] * v[4266] + v[3181] * v[4267] + v[3079] * v[4268]);
		v[4618] = v[3984] * v[4333];
		v[5288] = v[2217] * v[4332] + v[4618];
		v[5287] = v[2215] * v[4333] + v[4611];
		v[4334] = v[17] * (v[3293] * v[4266] + v[3183] * v[4267] + v[3081] * v[4268]);
		v[4623] = v[3982] * v[4334];
		v[5289] = v[2216] * v[4333] + v[4623];
		v[4620] = v[4334] * v[4755];
		v[4614] = v[4334] * v[4757];
		v[4335] = v[17] * (v[3295] * v[4266] + v[3185] * v[4267] + v[3083] * v[4268]);
		v[4643] = v[2212] * v[4335];
		v[4629] = v[3980] * v[4335];
		v[4336] = v[17] * (v[3297] * v[4266] + v[3187] * v[4267] + v[3085] * v[4268]);
		v[4636] = v[3978] * v[4336];
		v[5291] = v[2211] * v[4335] + v[4636];
		v[5290] = v[2209] * v[4336] + v[4629];
		v[4337] = v[17] * (v[3299] * v[4266] + v[3189] * v[4267] + v[3087] * v[4268]);
		v[4641] = v[3976] * v[4337];
		v[10756] = v[1409] * v[5238] + v[239] * v[5321];
		v[10757] = v[1409] * v[5237] + v[239] * v[5322];
		v[10758] = v[239] * v[5323];
		v[10759] = v[4641] + v[4336] * v[4750] + v[4335] * v[4752];
		v[10760] = v[2209] * v[4335] + v[2210] * v[4337] + v[4636];
		v[10761] = v[2211] * v[4336] + v[2212] * v[4337] + v[4629];
		v[10762] = v[1410] * v[5238] + v[240] * v[5321];
		v[10763] = v[1410] * v[5237] + v[240] * v[5322];
		v[10764] = v[240] * v[5323];
		v[10765] = v[4623] + v[4333] * v[4755] + v[4332] * v[4757];
		v[10766] = v[2215] * v[4332] + v[2216] * v[4334] + v[4618];
		v[10767] = v[2217] * v[4333] + v[2218] * v[4334] + v[4611];
		v[10768] = v[1411] * v[5238] - v[326] * v[5321];
		v[10769] = v[1411] * v[5237] - v[326] * v[5322];
		v[10770] = -(v[326] * v[5323]);
		v[10771] = v[4560] + v[4330] * v[4760] + v[4329] * v[4762];
		v[10772] = v[2221] * v[4329] + v[2222] * v[4331] + v[4555];
		v[10773] = v[2223] * v[4330] + v[2224] * v[4331] + v[4548];
		v[10774] = v[1412] * v[5238] - v[327] * v[5321];
		v[10775] = v[1412] * v[5237] - v[327] * v[5322];
		v[10776] = -(v[327] * v[5323]);
		v[10777] = v[4542] + v[4327] * v[4765] + v[4326] * v[4767];
		v[10778] = v[2227] * v[4326] + v[2228] * v[4328] + v[4537];
		v[10779] = v[2229] * v[4327] + v[2230] * v[4328] + v[4530];
		v[5292] = v[2210] * v[4336] + v[4641];
		v[4638] = v[4337] * v[4750];
		v[4632] = v[4337] * v[4752];
		v[4338] = v[4270] * v[985] + v[4273] * v[986];
		v[5312] = v[4338] + v[4278] * v[981];
		v[4339] = v[4274] + v[4276];
		v[4340] = v[4274] + v[4278];
		v[4341] = v[4276] + v[4278];
		v[4342] = v[5312] * v[962];
		v[4343] = v[4280] * v[981] + v[4281] * v[986];
		v[5313] = v[4343] + v[4274] * v[985];
		v[4344] = v[5313] * v[962];
		v[4345] = v[4283] * v[981] + v[4284] * v[985];
		v[5314] = v[4345] + v[4276] * v[986];
		v[4346] = v[4273] - v[4283] - v[4343] / 2e0 + v[4284] * v[982] + v[4270] * v[983] + v[4341] * v[984];
		v[4347] = -v[4281] + v[4284] - v[4338] / 2e0 + v[4283] * v[982] + v[4339] * v[983] + v[4280] * v[984];
		v[4348] = v[5314] * v[962];
		v[4349] = -v[4270] + v[4280] - v[4345] / 2e0 + v[4340] * v[982] + v[4273] * v[983] + v[4281] * v[984];
		v[4350] = v[4286] * v[979] + v[4288] * v[980];
		v[5308] = v[4350] + v[4292] * v[975];
		v[4351] = v[4289] + v[4290];
		v[4352] = v[4289] + v[4292];
		v[4353] = v[4290] + v[4292];
		v[4354] = v[5308] * v[931];
		v[4355] = v[4293] * v[975] + v[4294] * v[980];
		v[5309] = v[4355] + v[4289] * v[979];
		v[4356] = v[5309] * v[931];
		v[4357] = v[4295] * v[975] + v[4296] * v[979];
		v[5310] = v[4357] + v[4290] * v[980];
		v[4358] = v[4288] - v[4295] - v[4355] / 2e0 + v[4296] * v[976] + v[4286] * v[977] + v[4353] * v[978];
		v[4359] = -v[4294] + v[4296] - v[4350] / 2e0 + v[4295] * v[976] + v[4351] * v[977] + v[4293] * v[978];
		v[4360] = v[5310] * v[931];
		v[4361] = -v[4286] + v[4293] - v[4357] / 2e0 + v[4352] * v[976] + v[4288] * v[977] + v[4294] * v[978];
		v[4362] = v[4298] * v[973] + v[4300] * v[974];
		v[5304] = v[4362] + v[4304] * v[969];
		v[4363] = v[4301] + v[4302];
		v[4364] = v[4301] + v[4304];
		v[4365] = v[4302] + v[4304];
		v[4366] = v[5304] * v[780];
		v[4367] = v[4305] * v[969] + v[4306] * v[974];
		v[5305] = v[4367] + v[4301] * v[973];
		v[4368] = v[5305] * v[780];
		v[4369] = v[4307] * v[969] + v[4308] * v[973];
		v[5306] = v[4369] + v[4302] * v[974];
		v[4370] = v[4300] - v[4307] - v[4367] / 2e0 + v[4308] * v[970] + v[4298] * v[971] + v[4365] * v[972];
		v[4371] = -v[4306] + v[4308] - v[4362] / 2e0 + v[4307] * v[970] + v[4363] * v[971] + v[4305] * v[972];
		v[4372] = v[5306] * v[780];
		v[4373] = -v[4298] + v[4305] - v[4369] / 2e0 + v[4364] * v[970] + v[4300] * v[971] + v[4306] * v[972];
		v[4374] = v[4310] * v[967] + v[4312] * v[968];
		v[5300] = v[4374] + v[4316] * v[963];
		v[4375] = v[4313] + v[4314];
		v[4376] = v[4313] + v[4316];
		v[4377] = v[4314] + v[4316];
		v[4378] = v[5300] * v[749];
		v[4379] = v[4317] * v[963] + v[4318] * v[968];
		v[5301] = v[4379] + v[4313] * v[967];
		v[4380] = v[5301] * v[749];
		v[4381] = v[4319] * v[963] + v[4320] * v[967];
		v[5302] = v[4381] + v[4314] * v[968];
		v[4382] = v[4312] - v[4319] - v[4379] / 2e0 + v[4320] * v[964] + v[4310] * v[965] + v[4377] * v[966];
		v[4383] = -v[4318] + v[4320] - v[4374] / 2e0 + v[4319] * v[964] + v[4375] * v[965] + v[4317] * v[966];
		v[4384] = v[5302] * v[749];
		v[4385] = -v[4310] + v[4317] - v[4381] / 2e0 + v[4376] * v[964] + v[4312] * v[965] + v[4318] * v[966];
		v[4390] = v[17] * (v[4266] * v[4386] + v[4267] * v[4387]) + v[4268] * v[4389] + 2e0 * v[4388] * v[4720]
			+ v[4321] * v[5046];
		v[4396] = v[4267] * v[4394] + v[17] * (v[4266] * v[4391] + v[4268] * v[4395]) + 2e0 * v[4392] * v[4847]
			+ v[480] * v[5244];
		v[4402] = v[4266] * v[4399] + v[17] * (v[4267] * v[4400] + v[4268] * v[4401]) + 2e0 * v[4397] * v[4845]
			+ v[481] * v[5244];
		v[4403] = v[1483] * v[4220] + v[1484] * v[4223] + v[1485] * v[4226] + v[175] * v[4349] + v[4348] * v[4801]
			+ v[4342] * v[936] + v[4344] * v[944];
		v[4404] = v[1483] * v[4221] + v[1484] * v[4224] + v[1485] * v[4227] + v[175] * v[4346] + v[4344] * v[4800]
			+ v[4342] * v[935] + v[4348] * v[950];
		v[4405] = v[1483] * v[4222] + v[1484] * v[4225] + v[1485] * v[4228] + v[175] * v[4347] + v[4342] * v[4799]
			+ v[4344] * v[937] + v[4348] * v[945];
		v[4406] = v[1483] * v[4229] + v[1484] * v[4232] + v[1485] * v[4235] + v[194] * v[4361] + v[4360] * v[4798]
			+ v[4354] * v[905] + v[4356] * v[913];
		v[4407] = v[1483] * v[4230] + v[1484] * v[4233] + v[1485] * v[4236] + v[194] * v[4358] + v[4356] * v[4797]
			+ v[4354] * v[904] + v[4360] * v[919];
		v[4408] = v[1483] * v[4231] + v[1484] * v[4234] + v[1485] * v[4237] + v[194] * v[4359] + v[4354] * v[4796]
			+ v[4356] * v[906] + v[4360] * v[914];
		v[4409] = v[1483] * v[4242] + v[1484] * v[4245] + v[1485] * v[4248] + v[262] * v[4373] + v[4372] * v[4789]
			+ v[4366] * v[754] + v[4368] * v[762];
		v[4410] = v[1483] * v[4243] + v[1484] * v[4246] + v[1485] * v[4249] + v[262] * v[4370] + v[4368] * v[4788]
			+ v[4366] * v[753] + v[4372] * v[768];
		v[4411] = v[1483] * v[4244] + v[1484] * v[4247] + v[1485] * v[4250] + v[262] * v[4371] + v[4366] * v[4787]
			+ v[4368] * v[755] + v[4372] * v[763];
		v[4412] = v[1483] * v[4251] + v[1484] * v[4254] + v[1485] * v[4257] + v[281] * v[4385] + v[4384] * v[4786]
			+ v[4378] * v[723] + v[4380] * v[731];
		v[4413] = v[1483] * v[4252] + v[1484] * v[4255] + v[1485] * v[4258] + v[281] * v[4382] + v[4380] * v[4785]
			+ v[4378] * v[722] + v[4384] * v[737];
		v[4414] = v[1483] * v[4253] + v[1484] * v[4256] + v[1485] * v[4259] + v[281] * v[4383] + v[4378] * v[4784]
			+ v[4380] * v[724] + v[4384] * v[732];
		v[5246] = ((v[4402] * v[465] + v[4396] * v[466] + v[4390] * v[467]) * v[478] + (v[4268] * v[482] + v[5239]) * v[5245])
			/ v[1475];
		v[4417] = v[4390] * v[479] + v[467] * v[5246];
		v[4418] = v[4396] * v[479] + v[466] * v[5246];
		v[4419] = v[4402] * v[479] + v[465] * v[5246];
		v[4420] = -(v[326] * v[4417]) + v[4421];
		v[4422] = -(v[327] * v[4417]) - v[4421];
		v[4423] = v[239] * v[4417] + v[4424];
		v[4425] = v[240] * v[4417] - v[4424];
		v[4426] = -(v[326] * v[4418]) + v[4427];
		v[4428] = -(v[327] * v[4418]) - v[4427];
		v[4429] = v[239] * v[4418] + v[4430];
		v[4431] = v[240] * v[4418] - v[4430];
		v[4432] = -(v[326] * v[4419]) + v[4433];
		v[4434] = -(v[327] * v[4419]) - v[4433];
		v[4435] = v[239] * v[4419] + v[4436];
		v[4437] = v[240] * v[4419] - v[4436];
		v[4438] = v[2821] * v[4326];
		v[4439] = v[2820] * v[4326];
		v[4440] = v[2819] * v[4326];
		v[4441] = v[2816] * v[4327];
		v[4442] = v[4326] * v[451] + v[4327] * v[453];
		v[4443] = v[2814] * v[4327];
		v[4444] = -v[4438] + v[4441];
		v[4445] = v[4438] + v[4441];
		v[4446] = v[2815] * v[4327] + v[4442] * v[5247];
		v[4447] = v[2809] * v[4328];
		v[4448] = v[4326] * v[445] + v[4328] * v[453];
		v[4449] = v[4327] * v[445] + v[4328] * v[451];
		v[4450] = v[2811] * v[4328] + v[4448] * v[5248];
		v[4451] = v[4446] + v[4450];
		v[4452] = -v[4446] + v[4450];
		v[4453] = v[2810] * v[4328] + v[4449] * v[5249];
		v[4454] = v[4439] + v[4453];
		v[4455] = -v[4439] + v[4453];
		v[4456] = v[2806] * v[4329];
		v[4457] = v[2805] * v[4329];
		v[4458] = v[2804] * v[4329];
		v[4459] = v[2801] * v[4330];
		v[4460] = v[425] * v[4329] + v[427] * v[4330];
		v[4461] = v[2799] * v[4330];
		v[4462] = -v[4456] + v[4459];
		v[4463] = v[4456] + v[4459];
		v[4464] = v[2800] * v[4330] + v[4460] * v[5250];
		v[4465] = v[2794] * v[4331];
		v[4466] = v[419] * v[4329] + v[427] * v[4331];
		v[4467] = v[419] * v[4330] + v[425] * v[4331];
		v[4468] = v[2796] * v[4331] + v[4466] * v[5251];
		v[4469] = v[4464] + v[4468];
		v[4470] = -v[4464] + v[4468];
		v[4471] = v[2795] * v[4331] + v[4467] * v[5252];
		v[4472] = v[4457] + v[4471];
		v[4473] = -v[4457] + v[4471];
		v[4474] = v[2791] * v[4332];
		v[4475] = v[2790] * v[4332];
		v[4476] = v[2789] * v[4332];
		v[4477] = v[2786] * v[4333];
		v[4478] = v[399] * v[4332] + v[401] * v[4333];
		v[4479] = v[2784] * v[4333];
		v[4480] = -v[4474] + v[4477];
		v[4481] = v[4474] + v[4477];
		v[4482] = v[2785] * v[4333] + v[4478] * v[5253];
		v[4483] = v[2779] * v[4334];
		v[4484] = v[393] * v[4332] + v[401] * v[4334];
		v[4485] = v[393] * v[4333] + v[399] * v[4334];
		v[4486] = v[2781] * v[4334] + v[4484] * v[5254];
		v[4487] = v[4482] + v[4486];
		v[4488] = -v[4482] + v[4486];
		v[4489] = v[2780] * v[4334] + v[4485] * v[5255];
		v[4490] = v[4475] + v[4489];
		v[4491] = -v[4475] + v[4489];
		v[4492] = v[2776] * v[4335];
		v[4493] = v[2775] * v[4335];
		v[4494] = v[2774] * v[4335];
		v[4495] = v[2771] * v[4336];
		v[4496] = v[373] * v[4335] + v[375] * v[4336];
		v[4497] = v[2769] * v[4336];
		v[4498] = -v[4492] + v[4495];
		v[4499] = v[4492] + v[4495];
		v[4500] = v[2770] * v[4336] + v[4496] * v[5256];
		v[4501] = v[2764] * v[4337];
		v[4502] = v[367] * v[4335] + v[375] * v[4337];
		v[4503] = v[367] * v[4336] + v[373] * v[4337];
		v[4504] = v[2766] * v[4337] + v[4502] * v[5257];
		v[4505] = v[4500] + v[4504];
		v[4506] = -v[4500] + v[4504];
		v[4507] = v[2765] * v[4337] + v[4503] * v[5258];
		v[4508] = v[4493] + v[4507];
		v[4509] = -v[4493] + v[4507];
		v[4510] = v[339] * v[4432] + v[4263] * v[5259];
		v[4511] = v[337] * v[4432] + v[4261] * v[5259];
		v[4512] = v[339] * v[4434] + v[4263] * v[5260];
		v[4513] = v[337] * v[4434] + v[4261] * v[5260];
		v[4514] = v[339] * v[4426] + v[4263] * v[5261];
		v[4515] = v[337] * v[4426] + v[4261] * v[5261];
		v[4516] = v[339] * v[4428] + v[4263] * v[5262];
		v[4517] = v[337] * v[4428] + v[4261] * v[5262];
		v[4518] = v[339] * v[4420] + v[4263] * v[5263];
		v[4519] = v[337] * v[4420] + v[4261] * v[5263];
		v[4520] = v[326] * (-(v[149] * v[4409]) - v[152] * v[4410] - v[155] * v[4411]) + v[327] * (-(v[158] * v[4412])
			- v[161] * v[4413] - v[164] * v[4414]) + v[307] * v[4420] + v[316] * v[4422] + v[304] * v[4426] + v[313] * v[4428]
			+ v[301] * v[4432] + v[310] * v[4434];
		v[4521] = v[326] * (-(v[148] * v[4409]) - v[151] * v[4410] - v[154] * v[4411]) + v[327] * (-(v[157] * v[4412])
			- v[160] * v[4413] - v[163] * v[4414]) + v[306] * v[4420] + v[315] * v[4422] + v[303] * v[4426] + v[312] * v[4428]
			+ v[300] * v[4432] + v[309] * v[4434];
		v[4522] = v[339] * v[4422] + v[4263] * v[5264];
		v[4523] = v[337] * v[4422] + v[4261] * v[5264];
		v[5315] = v[328] * (v[4056] * v[4182] * v[5265] - v[5088] * (-(v[4520] * v[4748]) - v[4521] * v[5087]
			- v[4182] * v[5464]));
		v[4528] = (v[10059 + i4002] - v[10083 + i4002] + v[356] * v[4279] - v[359] * v[4279] + v[355] * v[4282] - v[358] * v[4282]
			+ v[3626] * v[4409] + v[3627] * v[4410] + v[3628] * v[4411] - v[3629] * v[4412] - v[3630] * v[4413] - v[3631] * v[4414]
			+ v[4324] * v[5090] - v[4325] * v[5090] + v[4419] * v[520] - v[4419] * v[521] + v[4418] * v[523] - v[4418] * v[524]
			+ v[4417] * v[526] + v[357] * v[5266] - v[360] * v[5266] - v[4526] * v[5267] + v[4527] * v[5267] - v[4417] * v[527]
			- v[4261] * v[5466] - v[4263] * v[5467] + v[4261] * v[5469] + v[4263] * v[5470]) / 2e0;
		v[4529] = v[164] * v[4522] + v[163] * v[4523] + v[453] * (v[4533] + v[5268]);
		v[4532] = v[161] * v[4522] + v[160] * v[4523] + v[4449] * v[4767] + v[451] * v[5268];
		v[4534] = v[158] * v[4522] + v[157] * v[4523] + v[445] * (v[4530] + v[4533]);
		v[4535] = v[164] * v[4516] + v[163] * v[4517] + v[4448] * v[4765] + v[453] * v[5269];
		v[4538] = v[161] * v[4516] + v[160] * v[4517] + v[451] * (v[4539] + v[5269]);
		v[4540] = v[158] * v[4516] + v[157] * v[4517] + v[445] * (v[4537] + v[4539]);
		v[4541] = v[2228] * v[4442] + v[164] * v[4512] + v[163] * v[4513] + v[453] * (v[4542] + v[4544]);
		v[4543] = v[161] * v[4512] + v[160] * v[4513] + v[451] * v[5270];
		v[4546] = v[158] * v[4512] + v[157] * v[4513] + v[445] * (v[4544] + v[5270]);
		v[4547] = v[155] * v[4518] + v[154] * v[4519] + v[427] * (v[4551] + v[5271]);
		v[4550] = v[152] * v[4518] + v[151] * v[4519] + v[4467] * v[4762] + v[425] * v[5271];
		v[4552] = v[149] * v[4518] + v[148] * v[4519] + v[419] * (v[4548] + v[4551]);
		v[4553] = v[155] * v[4514] + v[154] * v[4515] + v[4466] * v[4760] + v[427] * v[5272];
		v[4556] = v[152] * v[4514] + v[151] * v[4515] + v[425] * (v[4557] + v[5272]);
		v[4558] = v[149] * v[4514] + v[148] * v[4515] + v[419] * (v[4555] + v[4557]);
		v[4559] = v[2222] * v[4460] + v[155] * v[4510] + v[154] * v[4511] + v[427] * (v[4560] + v[4562]);
		v[4561] = v[152] * v[4510] + v[151] * v[4511] + v[425] * v[5273];
		v[4564] = v[149] * v[4510] + v[148] * v[4511] + v[419] * (v[4562] + v[5273]);
		v[4565] = -(v[281] * v[4444]) - v[298] * v[4451] + v[297] * v[4454] + v[4447] * v[4993];
		v[4566] = -(v[297] * v[4445]) - v[299] * v[4451] - v[281] * v[4455] + v[4443] * v[4992];
		v[4567] = -(v[298] * v[4445]) - v[281] * v[4452] + v[299] * v[4454] + v[4440] * v[4991];
		v[4568] = v[4004] * v[4572] + v[4529] * v[4741] + v[4378] * v[5274];
		v[4569] = v[1348] * v[4378] + v[281] * v[4532] + v[4005] * v[5236];
		v[4570] = v[1349] * v[4378] + v[281] * v[4534] + v[4006] * v[5236];
		v[4571] = v[1347] * v[4380] + v[281] * v[4535] + v[4016] * v[5236];
		v[4574] = v[1349] * v[4380] + v[281] * v[4540] + v[4018] * v[5236];
		v[4575] = v[1347] * v[4384] + v[281] * v[4541] + v[4036] * v[5236];
		v[4576] = v[1348] * v[4384] + v[281] * v[4543] + v[4037] * v[5236];
		v[4577] = v[4038] * v[4572] + v[4546] * v[4741] + v[4384] * v[5276];
		v[4578] = -(v[262] * v[4462]) - v[279] * v[4469] + v[278] * v[4472] + v[4465] * v[4989];
		v[4579] = -(v[278] * v[4463]) - v[280] * v[4469] - v[262] * v[4473] + v[4461] * v[4988];
		v[4580] = -(v[279] * v[4463]) - v[262] * v[4470] + v[280] * v[4472] + v[4458] * v[4987];
		v[4581] = v[4007] * v[4585] + v[4547] * v[4737] + v[4366] * v[5277];
		v[4582] = v[1351] * v[4366] + v[262] * v[4550] + v[4008] * v[5235];
		v[4583] = v[1352] * v[4366] + v[262] * v[4552] + v[4009] * v[5235];
		v[4584] = v[1350] * v[4368] + v[262] * v[4553] + v[4019] * v[5235];
		v[4587] = v[1352] * v[4368] + v[262] * v[4558] + v[4021] * v[5235];
		v[4588] = v[1350] * v[4372] + v[262] * v[4559] + v[4039] * v[5235];
		v[4589] = v[1351] * v[4372] + v[262] * v[4561] + v[4040] * v[5235];
		v[4590] = v[4041] * v[4585] + v[4564] * v[4737] + v[4372] * v[5279];
		v[4591] = v[252] * v[4435] + v[4241] * v[5280];
		v[4592] = v[250] * v[4435] + v[4239] * v[5280];
		v[4593] = v[252] * v[4437] + v[4241] * v[5281];
		v[4594] = v[250] * v[4437] + v[4239] * v[5281];
		v[4595] = v[252] * v[4429] + v[4241] * v[5282];
		v[4596] = v[250] * v[4429] + v[4239] * v[5282];
		v[4597] = v[252] * v[4431] + v[4241] * v[5283];
		v[4598] = v[250] * v[4431] + v[4239] * v[5283];
		v[4599] = v[252] * v[4423] + v[4241] * v[5284];
		v[4600] = v[250] * v[4423] + v[4239] * v[5284];
		v[4601] = v[220] * v[4423] + v[229] * v[4425] + v[217] * v[4429] + v[226] * v[4431] + v[214] * v[4435] + v[223] * v[4437]
			+ v[239] * (v[4403] * v[86] + v[4404] * v[89] + v[4405] * v[92]) + v[240] * (v[101] * v[4408] + v[4406] * v[95]
				+ v[4407] * v[98]);
		v[4602] = v[219] * v[4423] + v[228] * v[4425] + v[216] * v[4429] + v[225] * v[4431] + v[213] * v[4435] + v[222] * v[4437]
			+ v[239] * (v[4403] * v[85] + v[4404] * v[88] + v[4405] * v[91]) + v[240] * (v[100] * v[4408] + v[4406] * v[94]
				+ v[4407] * v[97]);
		v[4603] = v[252] * v[4425] + v[4241] * v[5285];
		v[4604] = v[250] * v[4425] + v[4239] * v[5285];
		v[5316] = v[241] * (v[4078] * v[4164] * v[5286] - v[5106] * (-((-(v[1457] * v[4049] * v[4164]) + v[1564] * v[4602]
			) * v[76]) - (v[1564] * v[4045] * v[4164] + v[1457] * v[4601]) * v[77]));
		v[4609] = (v[10107 + i4002] - v[10131 + i4002] - v[350] * v[4279] + v[353] * v[4279] - v[349] * v[4282] + v[352] * v[4282]
			- v[3760] * v[4403] - v[3761] * v[4404] - v[3762] * v[4405] + v[3763] * v[4406] + v[3764] * v[4407] + v[3765] * v[4408]
			- v[4419] * v[492] + v[4419] * v[493] - v[4418] * v[495] + v[4418] * v[496] - v[4417] * v[498] + v[4417] * v[499]
			- v[4322] * v[5090] + v[4323] * v[5090] - v[351] * v[5266] + v[354] * v[5266] - v[4607] * v[5267] + v[4608] * v[5267]
			- v[4239] * v[5472] - v[4241] * v[5473] + v[4239] * v[5475] + v[4241] * v[5476]) / 2e0;
		v[4610] = v[101] * v[4603] + v[100] * v[4604] + v[401] * (v[4614] + v[5287]);
		v[4613] = v[4485] * v[4757] + v[399] * v[5287] + v[4604] * v[97] + v[4603] * v[98];
		v[4615] = v[393] * (v[4611] + v[4614]) + v[4604] * v[94] + v[4603] * v[95];
		v[4616] = v[101] * v[4597] + v[100] * v[4598] + v[4484] * v[4755] + v[401] * v[5288];
		v[4619] = v[399] * (v[4620] + v[5288]) + v[4598] * v[97] + v[4597] * v[98];
		v[4621] = v[393] * (v[4618] + v[4620]) + v[4598] * v[94] + v[4597] * v[95];
		v[4622] = v[2216] * v[4478] + v[101] * v[4593] + v[100] * v[4594] + v[401] * (v[4623] + v[4625]);
		v[4624] = v[399] * v[5289] + v[4594] * v[97] + v[4593] * v[98];
		v[4627] = v[393] * (v[4625] + v[5289]) + v[4594] * v[94] + v[4593] * v[95];
		v[4628] = v[375] * (v[4632] + v[5290]) + v[4600] * v[91] + v[4599] * v[92];
		v[4631] = v[4503] * v[4752] + v[373] * v[5290] + v[4600] * v[88] + v[4599] * v[89];
		v[4633] = v[367] * (v[4629] + v[4632]) + v[4600] * v[85] + v[4599] * v[86];
		v[4634] = v[4502] * v[4750] + v[375] * v[5291] + v[4596] * v[91] + v[4595] * v[92];
		v[4637] = v[373] * (v[4638] + v[5291]) + v[4596] * v[88] + v[4595] * v[89];
		v[4639] = v[367] * (v[4636] + v[4638]) + v[4596] * v[85] + v[4595] * v[86];
		v[4640] = v[2210] * v[4496] + v[375] * (v[4641] + v[4643]) + v[4592] * v[91] + v[4591] * v[92];
		v[4642] = v[373] * v[5292] + v[4592] * v[88] + v[4591] * v[89];
		v[4645] = v[367] * (v[4643] + v[5292]) + v[4592] * v[85] + v[4591] * v[86];
		v[4646] = -(v[194] * v[4480]) - v[211] * v[4487] + v[210] * v[4490] + v[4483] * v[4979];
		v[4647] = -(v[210] * v[4481]) - v[212] * v[4487] - v[194] * v[4491] + v[4479] * v[4978];
		v[4648] = -(v[211] * v[4481]) - v[194] * v[4488] + v[212] * v[4490] + v[4476] * v[4977];
		v[4649] = v[4010] * v[4653] + v[4610] * v[4728] + v[4354] * v[5293];
		v[4650] = v[1354] * v[4354] + v[194] * v[4613] + v[4011] * v[5234];
		v[4651] = v[1355] * v[4354] + v[194] * v[4615] + v[4012] * v[5234];
		v[4652] = v[1353] * v[4356] + v[194] * v[4616] + v[4022] * v[5234];
		v[4655] = v[1355] * v[4356] + v[194] * v[4621] + v[4024] * v[5234];
		v[4656] = v[1353] * v[4360] + v[194] * v[4622] + v[4050] * v[5234];
		v[4657] = v[1354] * v[4360] + v[194] * v[4624] + v[4051] * v[5234];
		v[4658] = v[4052] * v[4653] + v[4627] * v[4728] + v[4360] * v[5295];
		v[4659] = -(v[175] * v[4498]) - v[192] * v[4505] + v[191] * v[4508] + v[4501] * v[4975];
		v[4660] = -(v[191] * v[4499]) - v[193] * v[4505] - v[175] * v[4509] + v[4497] * v[4974];
		v[4661] = -(v[192] * v[4499]) - v[175] * v[4506] + v[193] * v[4508] + v[4494] * v[4973];
		v[4662] = v[4013] * v[4666] + v[4628] * v[4724] + v[4342] * v[5296];
		v[4663] = v[1357] * v[4342] + v[175] * v[4631] + v[4014] * v[5233];
		v[4664] = v[1358] * v[4342] + v[175] * v[4633] + v[4015] * v[5233];
		v[4665] = v[1356] * v[4344] + v[175] * v[4634] + v[4025] * v[5233];
		v[4668] = v[1358] * v[4344] + v[175] * v[4639] + v[4027] * v[5233];
		v[4669] = v[1356] * v[4348] + v[175] * v[4640] + v[4053] * v[5233];
		v[4670] = v[1357] * v[4348] + v[175] * v[4642] + v[4054] * v[5233];
		v[4671] = v[4055] * v[4666] + v[4645] * v[4724] + v[4348] * v[5298];
		v[4674] = v[5299] * (v[4175] * v[4673] + v[1364] * v[5300] + v[1363] * v[5301] + v[1362] * v[5302]) + v[749] * (
			(v[4038] * v[4177] + v[4017] * v[4179] + v[4004] * v[4181]) / 2e0 + v[4037] * v[4211] + v[4018] * v[4213]
			+ v[4016] * v[4214] + v[4005] * v[4216] + v[4036] * v[4217] + v[4006] * v[4219] + v[1348] * v[4382] + v[1347] * v[4383]
			+ v[1349] * v[4385] - v[299] * v[4444] - v[297] * v[4452] - v[298] * v[4455] + v[4529] * v[4784] + v[4538] * v[4785]
			+ v[4546] * v[4786] + (v[2818] * v[4326] + v[2813] * v[4327] + v[2808] * v[4328] + v[4440] + v[2616] * v[4442] + v[4443]
				+ v[4447] + v[2614] * v[4448] + v[2610] * v[4449]) * v[5113] + v[4532] * v[722] + v[4534] * v[723] + v[4535] * v[724]
			+ v[4540] * v[731] + v[4541] * v[732] + v[4543] * v[737] - v[4565] * v[964] - v[4567] * v[965] + v[4566] * v[966]);
		v[5320] = -(v[4017] * v[4572]) + v[4674] - v[4538] * v[4741] - v[4380] * v[5275];
		v[4675] = v[4570] + v[4575];
		v[4676] = v[4569] + v[4571];
		v[4677] = v[4574] + v[4576];
		v[4680] = v[5303] * (v[4167] * v[4679] + v[1370] * v[5304] + v[1369] * v[5305] + v[1368] * v[5306]) + v[780] * (
			(v[4041] * v[4169] + v[4020] * v[4171] + v[4007] * v[4173]) / 2e0 + v[4040] * v[4202] + v[4021] * v[4204]
			+ v[4019] * v[4205] + v[4008] * v[4207] + v[4039] * v[4208] + v[4009] * v[4210] + v[1351] * v[4370] + v[1350] * v[4371]
			+ v[1352] * v[4373] - v[280] * v[4462] - v[278] * v[4470] - v[279] * v[4473] + v[4547] * v[4787] + v[4556] * v[4788]
			+ v[4564] * v[4789] + (v[2803] * v[4329] + v[2798] * v[4330] + v[2793] * v[4331] + v[4458] + v[2588] * v[4460] + v[4461]
				+ v[4465] + v[2586] * v[4466] + v[2582] * v[4467]) * v[5126] + v[4550] * v[753] + v[4552] * v[754] + v[4553] * v[755]
			+ v[4558] * v[762] + v[4559] * v[763] + v[4561] * v[768] - v[4578] * v[970] - v[4580] * v[971] + v[4579] * v[972]);
		v[5319] = -(v[4020] * v[4585]) + v[4680] - v[4556] * v[4737] - v[4368] * v[5278];
		v[4681] = v[4583] + v[4588];
		v[4682] = v[4582] + v[4584];
		v[4683] = v[4587] + v[4589];
		v[4686] = v[5307] * (v[4157] * v[4685] + v[1376] * v[5308] + v[1375] * v[5309] + v[1374] * v[5310]) + v[931] * (
			(v[4052] * v[4159] + v[4023] * v[4161] + v[4010] * v[4163]) / 2e0 + v[4051] * v[4193] + v[4024] * v[4195]
			+ v[4022] * v[4196] + v[4011] * v[4198] + v[4050] * v[4199] + v[4012] * v[4201] + v[1354] * v[4358] + v[1353] * v[4359]
			+ v[1355] * v[4361] - v[212] * v[4480] - v[210] * v[4488] - v[211] * v[4491] + v[4610] * v[4796] + v[4619] * v[4797]
			+ v[4627] * v[4798] + (v[2788] * v[4332] + v[2783] * v[4333] + v[2778] * v[4334] + v[4476] + v[2560] * v[4478] + v[4479]
				+ v[4483] + v[2558] * v[4484] + v[2554] * v[4485]) * v[5139] + v[4613] * v[904] + v[4615] * v[905] + v[4616] * v[906]
			+ v[4621] * v[913] + v[4622] * v[914] + v[4624] * v[919] - v[4646] * v[976] - v[4648] * v[977] + v[4647] * v[978]);
		v[5318] = -(v[4023] * v[4653]) + v[4686] - v[4619] * v[4728] - v[4356] * v[5294];
		v[4687] = v[4651] + v[4656];
		v[4688] = v[4650] + v[4652];
		v[4689] = v[4655] + v[4657];
		v[4692] = v[5311] * (v[4149] * v[4691] + v[1382] * v[5312] + v[1381] * v[5313] + v[1380] * v[5314]) + v[962] * (
			(v[4055] * v[4151] + v[4026] * v[4153] + v[4013] * v[4155]) / 2e0 + v[4054] * v[4184] + v[4027] * v[4186]
			+ v[4025] * v[4187] + v[4014] * v[4189] + v[4053] * v[4190] + v[4015] * v[4192] + v[1357] * v[4346] + v[1356] * v[4347]
			+ v[1358] * v[4349] - v[193] * v[4498] - v[191] * v[4506] - v[192] * v[4509] + v[4628] * v[4799] + v[4637] * v[4800]
			+ v[4645] * v[4801] + (v[2773] * v[4335] + v[2768] * v[4336] + v[2763] * v[4337] + v[4494] + v[2532] * v[4496] + v[4497]
				+ v[4501] + v[2530] * v[4502] + v[2526] * v[4503]) * v[5152] + v[4631] * v[935] + v[4633] * v[936] + v[4634] * v[937]
			+ v[4639] * v[944] + v[4640] * v[945] + v[4642] * v[950] - v[4659] * v[982] - v[4661] * v[983] + v[4660] * v[984]);
		v[5317] = -(v[4026] * v[4666]) + v[4692] - v[4637] * v[4724] - v[4344] * v[5297];
		v[4693] = v[4664] + v[4669];
		v[4694] = v[4663] + v[4665];
		v[4695] = v[4668] + v[4670];
		v[10684] = 0e0;
		v[10685] = 0e0;
		v[10686] = 0e0;
		v[10687] = 0e0;
		v[10688] = v[4111];
		v[10689] = v[4109];
		v[10690] = 0e0;
		v[10691] = 0e0;
		v[10692] = 0e0;
		v[10693] = 0e0;
		v[10694] = 0e0;
		v[10695] = 0e0;
		v[10696] = 0e0;
		v[10697] = 0e0;
		v[10698] = 0e0;
		v[10699] = 0e0;
		v[10700] = 0e0;
		v[10701] = 0e0;
		v[10702] = 0e0;
		v[10703] = 0e0;
		v[10704] = 0e0;
		v[10705] = 0e0;
		v[10706] = 0e0;
		v[10707] = 0e0;
		v[10636] = 0e0;
		v[10637] = 0e0;
		v[10638] = 0e0;
		v[10639] = v[4111];
		v[10640] = 0e0;
		v[10641] = v[4110];
		v[10642] = 0e0;
		v[10643] = 0e0;
		v[10644] = 0e0;
		v[10645] = 0e0;
		v[10646] = 0e0;
		v[10647] = 0e0;
		v[10648] = 0e0;
		v[10649] = 0e0;
		v[10650] = 0e0;
		v[10651] = 0e0;
		v[10652] = 0e0;
		v[10653] = 0e0;
		v[10654] = 0e0;
		v[10655] = 0e0;
		v[10656] = 0e0;
		v[10657] = 0e0;
		v[10658] = 0e0;
		v[10659] = 0e0;
		v[10612] = 0e0;
		v[10613] = 0e0;
		v[10614] = 0e0;
		v[10615] = v[4109];
		v[10616] = v[4110];
		v[10617] = 0e0;
		v[10618] = 0e0;
		v[10619] = 0e0;
		v[10620] = 0e0;
		v[10621] = 0e0;
		v[10622] = 0e0;
		v[10623] = 0e0;
		v[10624] = 0e0;
		v[10625] = 0e0;
		v[10626] = 0e0;
		v[10627] = 0e0;
		v[10628] = 0e0;
		v[10629] = 0e0;
		v[10630] = 0e0;
		v[10631] = 0e0;
		v[10632] = 0e0;
		v[10633] = 0e0;
		v[10634] = 0e0;
		v[10635] = 0e0;
		v[10540] = 0e0;
		v[10541] = 0e0;
		v[10542] = 0e0;
		v[10543] = 0e0;
		v[10544] = 0e0;
		v[10545] = 0e0;
		v[10546] = 0e0;
		v[10547] = 0e0;
		v[10548] = 0e0;
		v[10549] = 0e0;
		v[10550] = v[4108];
		v[10551] = v[4106];
		v[10552] = 0e0;
		v[10553] = 0e0;
		v[10554] = 0e0;
		v[10555] = 0e0;
		v[10556] = 0e0;
		v[10557] = 0e0;
		v[10558] = 0e0;
		v[10559] = 0e0;
		v[10560] = 0e0;
		v[10561] = 0e0;
		v[10562] = 0e0;
		v[10563] = 0e0;
		v[10492] = 0e0;
		v[10493] = 0e0;
		v[10494] = 0e0;
		v[10495] = 0e0;
		v[10496] = 0e0;
		v[10497] = 0e0;
		v[10498] = 0e0;
		v[10499] = 0e0;
		v[10500] = 0e0;
		v[10501] = v[4108];
		v[10502] = 0e0;
		v[10503] = v[4107];
		v[10504] = 0e0;
		v[10505] = 0e0;
		v[10506] = 0e0;
		v[10507] = 0e0;
		v[10508] = 0e0;
		v[10509] = 0e0;
		v[10510] = 0e0;
		v[10511] = 0e0;
		v[10512] = 0e0;
		v[10513] = 0e0;
		v[10514] = 0e0;
		v[10515] = 0e0;
		v[10468] = 0e0;
		v[10469] = 0e0;
		v[10470] = 0e0;
		v[10471] = 0e0;
		v[10472] = 0e0;
		v[10473] = 0e0;
		v[10474] = 0e0;
		v[10475] = 0e0;
		v[10476] = 0e0;
		v[10477] = v[4106];
		v[10478] = v[4107];
		v[10479] = 0e0;
		v[10480] = 0e0;
		v[10481] = 0e0;
		v[10482] = 0e0;
		v[10483] = 0e0;
		v[10484] = 0e0;
		v[10485] = 0e0;
		v[10486] = 0e0;
		v[10487] = 0e0;
		v[10488] = 0e0;
		v[10489] = 0e0;
		v[10490] = 0e0;
		v[10491] = 0e0;
		v[10396] = 0e0;
		v[10397] = 0e0;
		v[10398] = 0e0;
		v[10399] = 0e0;
		v[10400] = 0e0;
		v[10401] = 0e0;
		v[10402] = 0e0;
		v[10403] = 0e0;
		v[10404] = 0e0;
		v[10405] = 0e0;
		v[10406] = 0e0;
		v[10407] = 0e0;
		v[10408] = 0e0;
		v[10409] = 0e0;
		v[10410] = 0e0;
		v[10411] = 0e0;
		v[10412] = v[4105];
		v[10413] = v[4103];
		v[10414] = 0e0;
		v[10415] = 0e0;
		v[10416] = 0e0;
		v[10417] = 0e0;
		v[10418] = 0e0;
		v[10419] = 0e0;
		v[10348] = 0e0;
		v[10349] = 0e0;
		v[10350] = 0e0;
		v[10351] = 0e0;
		v[10352] = 0e0;
		v[10353] = 0e0;
		v[10354] = 0e0;
		v[10355] = 0e0;
		v[10356] = 0e0;
		v[10357] = 0e0;
		v[10358] = 0e0;
		v[10359] = 0e0;
		v[10360] = 0e0;
		v[10361] = 0e0;
		v[10362] = 0e0;
		v[10363] = v[4105];
		v[10364] = 0e0;
		v[10365] = v[4104];
		v[10366] = 0e0;
		v[10367] = 0e0;
		v[10368] = 0e0;
		v[10369] = 0e0;
		v[10370] = 0e0;
		v[10371] = 0e0;
		v[10324] = 0e0;
		v[10325] = 0e0;
		v[10326] = 0e0;
		v[10327] = 0e0;
		v[10328] = 0e0;
		v[10329] = 0e0;
		v[10330] = 0e0;
		v[10331] = 0e0;
		v[10332] = 0e0;
		v[10333] = 0e0;
		v[10334] = 0e0;
		v[10335] = 0e0;
		v[10336] = 0e0;
		v[10337] = 0e0;
		v[10338] = 0e0;
		v[10339] = v[4103];
		v[10340] = v[4104];
		v[10341] = 0e0;
		v[10342] = 0e0;
		v[10343] = 0e0;
		v[10344] = 0e0;
		v[10345] = 0e0;
		v[10346] = 0e0;
		v[10347] = 0e0;
		v[10252] = 0e0;
		v[10253] = 0e0;
		v[10254] = 0e0;
		v[10255] = 0e0;
		v[10256] = 0e0;
		v[10257] = 0e0;
		v[10258] = 0e0;
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
		v[10269] = 0e0;
		v[10270] = 0e0;
		v[10271] = 0e0;
		v[10272] = 0e0;
		v[10273] = 0e0;
		v[10274] = v[4102];
		v[10275] = v[4100];
		v[10204] = 0e0;
		v[10205] = 0e0;
		v[10206] = 0e0;
		v[10207] = 0e0;
		v[10208] = 0e0;
		v[10209] = 0e0;
		v[10210] = 0e0;
		v[10211] = 0e0;
		v[10212] = 0e0;
		v[10213] = 0e0;
		v[10214] = 0e0;
		v[10215] = 0e0;
		v[10216] = 0e0;
		v[10217] = 0e0;
		v[10218] = 0e0;
		v[10219] = 0e0;
		v[10220] = 0e0;
		v[10221] = 0e0;
		v[10222] = 0e0;
		v[10223] = 0e0;
		v[10224] = 0e0;
		v[10225] = v[4102];
		v[10226] = 0e0;
		v[10227] = v[4101];
		v[10180] = 0e0;
		v[10181] = 0e0;
		v[10182] = 0e0;
		v[10183] = 0e0;
		v[10184] = 0e0;
		v[10185] = 0e0;
		v[10186] = 0e0;
		v[10187] = 0e0;
		v[10188] = 0e0;
		v[10189] = 0e0;
		v[10190] = 0e0;
		v[10191] = 0e0;
		v[10192] = 0e0;
		v[10193] = 0e0;
		v[10194] = 0e0;
		v[10195] = 0e0;
		v[10196] = 0e0;
		v[10197] = 0e0;
		v[10198] = 0e0;
		v[10199] = 0e0;
		v[10200] = 0e0;
		v[10201] = v[4100];
		v[10202] = v[4101];
		v[10203] = 0e0;
		v[10732] = v[4435];
		v[10733] = v[4429];
		v[10734] = v[4423];
		v[10735] = (v[10683 + i4002] + v[1377] * v[4283] + v[1378] * v[4284] + v[1379] * v[4340]) / 2e0 + v[4663] - v[4665]
			- v[4659] * v[4724] + 2e0 * (v[1407] * v[4273] + v[1403] * v[4281] + v[4188] * v[4708] + v[4276] * v[4843])
			+ v[4693] * v[983] + v[4695] * v[984] + (-v[4662] + v[5317]) * v[986];
		v[10736] = (v[10635 + i4002] + v[1377] * v[4280] + v[1379] * v[4281] + v[1378] * v[4341]) / 2e0 - v[4664] + v[4669]
			+ v[4660] * v[4724] + 2e0 * (v[1407] * v[4270] + v[1402] * v[4284] + v[4191] * v[4709] + v[4274] * v[4846])
			+ v[4695] * v[982] + v[4694] * v[983] + (-v[4662] - v[4671] + v[4692]) * v[985];
		v[10737] = (v[10611 + i4002] + v[1378] * v[4270] + v[1379] * v[4273] + v[1377] * v[4339]) / 2e0 + v[4668] - v[4670]
			- v[4661] * v[4724] + 2e0 * (v[1403] * v[4280] + v[1402] * v[4283] + v[4185] * v[4710] + v[4278] * v[4859]) + (-v[4671]
				+ v[5317]) * v[981] + v[4693] * v[982] + v[4694] * v[984];
		v[10738] = v[4437];
		v[10739] = v[4431];
		v[10740] = v[4425];
		v[10741] = (v[10539 + i4002] + v[1371] * v[4295] + v[1372] * v[4296] + v[1373] * v[4352]) / 2e0 + v[4650] - v[4652]
			- v[4646] * v[4728] + 2e0 * (v[1401] * v[4288] + v[1397] * v[4294] + v[4197] * v[4711] + v[4290] * v[4837])
			+ v[4687] * v[977] + v[4689] * v[978] + (-v[4649] + v[5318]) * v[980];
		v[10742] = (v[10491 + i4002] + v[1371] * v[4293] + v[1373] * v[4294] + v[1372] * v[4353]) / 2e0 - v[4651] + v[4656]
			+ v[4647] * v[4728] + 2e0 * (v[1401] * v[4286] + v[1396] * v[4296] + v[4200] * v[4712] + v[4289] * v[4838])
			+ v[4689] * v[976] + v[4688] * v[977] + (-v[4649] - v[4658] + v[4686]) * v[979];
		v[10743] = (v[10467 + i4002] + v[1372] * v[4286] + v[1373] * v[4288] + v[1371] * v[4351]) / 2e0 + v[4655] - v[4657]
			- v[4648] * v[4728] + 2e0 * (v[1397] * v[4293] + v[1396] * v[4295] + v[4194] * v[4713] + v[4292] * v[4839]) + (-v[4658]
				+ v[5318]) * v[975] + v[4687] * v[976] + v[4688] * v[978];
		v[10744] = v[4432];
		v[10745] = v[4426];
		v[10746] = v[4420];
		v[10747] = (v[10395 + i4002] + v[1365] * v[4307] + v[1366] * v[4308] + v[1367] * v[4364]) / 2e0 + v[4582] - v[4584]
			- v[4578] * v[4737] + 2e0 * (v[1395] * v[4300] + v[1391] * v[4306] + v[4206] * v[4714] + v[4302] * v[4831])
			+ v[4681] * v[971] + v[4683] * v[972] + (-v[4581] + v[5319]) * v[974];
		v[10748] = (v[10347 + i4002] + v[1365] * v[4305] + v[1367] * v[4306] + v[1366] * v[4365]) / 2e0 - v[4583] + v[4588]
			+ v[4579] * v[4737] + 2e0 * (v[1395] * v[4298] + v[1390] * v[4308] + v[4209] * v[4715] + v[4301] * v[4832])
			+ v[4683] * v[970] + v[4682] * v[971] + (-v[4581] - v[4590] + v[4680]) * v[973];
		v[10749] = (v[10323 + i4002] + v[1366] * v[4298] + v[1367] * v[4300] + v[1365] * v[4363]) / 2e0 + v[4587] - v[4589]
			- v[4580] * v[4737] + 2e0 * (v[1391] * v[4305] + v[1390] * v[4307] + v[4203] * v[4716] + v[4304] * v[4833]) + (-v[4590]
				+ v[5319]) * v[969] + v[4681] * v[970] + v[4682] * v[972];
		v[10750] = v[4434];
		v[10751] = v[4428];
		v[10752] = v[4422];
		v[10753] = (v[10251 + i4002] + v[1359] * v[4319] + v[1360] * v[4320] + v[1361] * v[4376]) / 2e0 + v[4569] - v[4571]
			- v[4565] * v[4741] + 2e0 * (v[1389] * v[4312] + v[1385] * v[4318] + v[4215] * v[4717] + v[4314] * v[4825])
			+ v[4675] * v[965] + v[4677] * v[966] + (-v[4568] + v[5320]) * v[968];
		v[10754] = (v[10203 + i4002] + v[1359] * v[4317] + v[1361] * v[4318] + v[1360] * v[4377]) / 2e0 - v[4570] + v[4575]
			+ v[4566] * v[4741] + 2e0 * (v[1389] * v[4310] + v[1384] * v[4320] + v[4218] * v[4718] + v[4313] * v[4826])
			+ v[4677] * v[964] + v[4676] * v[965] + (-v[4568] - v[4577] + v[4674]) * v[967];
		v[10755] = (v[10179 + i4002] + v[1360] * v[4310] + v[1361] * v[4312] + v[1359] * v[4375]) / 2e0 + v[4574] - v[4576]
			- v[4567] * v[4741] + 2e0 * (v[1385] * v[4317] + v[1384] * v[4319] + v[4212] * v[4719] + v[4316] * v[4827]) + (-v[4577]
				+ v[5320]) * v[963] + v[4675] * v[964] + v[4676] * v[966];
		v[4698] = -(v[4182] * v[5225]);
		v[4701] = v[1565] * (v[140] * (v[4031] * v[4183] + v[335] * v[4520]) - v[4182] * v[4700] + v[1463] * (v[4697] * v[4698]
			+ v[3966] * v[5315])) - v[1462] * (v[139] * (v[4035] * v[4183] + v[335] * v[4521]) + v[4182] * v[4696] + v[1465] *
				(v[4698] * v[4699] + v[3968] * v[5315]));
		v[4704] = -(v[4164] * v[5229]);
		v[4707] = -(v[1457] * (v[4164] * v[4702] + v[1460] * (v[4704] * v[4705] + v[3973] * v[5316]) + (v[4049] * v[4165]
			+ v[248] * v[4602]) * v[76])) + v[1564] * (-(v[4164] * v[4706]) + v[1458] * (v[4703] * v[4704] + v[3971] * v[5316]) +
				(v[4045] * v[4165] + v[248] * v[4601]) * v[77]);
		Rc[i4002 - 1] += v[4115] * v[4164] + v[4113] * v[4182] + v[4079] * v[4264] + v[4057] * v[4265] + v[9932 + i4002];
		for (i4122 = 1; i4122 <= 24; i4122++) {
			Kc[i4002 - 1][i4122 - 1] += v[10731 + i4122] + v[10755 + i4122] * v[19] + v[4707] * v[9829 + i4122] + v[4701] * v[9853
				+ i4122] + v[4609] * v[9877 + i4122] + v[4528] * v[9901 + i4122];
		};/* end for */
	};/* end for */
#pragma endregion

}