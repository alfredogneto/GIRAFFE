#include "SplineElement_SplineElement.h"

#include "SplineElement.h"
#include "Dynamic.h"
#include "Spline.h"
#include "SPContactData.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

SplineElement_SplineElement::SplineElement_SplineElement()
{
	DefaultValues();
}

SplineElement_SplineElement::~SplineElement_SplineElement()
{
	Free();
}

void SplineElement_SplineElement::InitializeConvectiveRange()
{
	SplineElement* surf1;		//Ponteiro para a superfície 1
	SplineElement* surf2;		//Ponteiro para a superfície 2
	surf1 = static_cast<SplineElement*>(db.splines[spline1_ID - 1]->sp_element[surf1_ID]);
	surf2 = static_cast<SplineElement*>(db.splines[spline2_ID - 1]->sp_element[surf2_ID]);

	convective_min(0, 0) = surf1->knot_element[2];
	convective_max(0, 0) = surf1->knot_element[3];
	convective_range(0, 0) = surf1->knot_element[3] - surf1->knot_element[2];

	convective_min(1, 0) = surf2->knot_element[2];
	convective_max(1, 0) = surf2->knot_element[3];
	convective_range(1, 0) = surf2->knot_element[3] - surf2->knot_element[2];

}

//Verifica range de coordenadas convectivas
int SplineElement_SplineElement::VerifyConvectiveRange(Matrix& mc)
{
	int return_value;

	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	SplineElement* surf1;		//Ponteiro para a superfície 1
	SplineElement* surf2;		//Ponteiro para a superfície 2
	surf1 = static_cast<SplineElement*>(db.splines[spline1_ID - 1]->sp_element[surf1_ID]);
	surf2 = static_cast<SplineElement*>(db.splines[spline2_ID - 1]->sp_element[surf2_ID]);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Intervalo de coordenadas convectivas definidas pelo knot vector e range de coordenadas convectivas para cada um dos trechos de spline
	double surf1_knotA = surf1->knot_element[2];
	double surf1_knotB = surf1->knot_element[3];
	double surf1_range = surf1_knotB - surf1_knotA;

	double surf2_knotA = surf2->knot_element[2];
	double surf2_knotB = surf2->knot_element[3];
	double surf2_range = surf2_knotB - surf2_knotA;

	double perc = 0.1;



	//Retornos:
	//0 - Range fisico da superficie
	//4 - Fora do range fisico da superficie - proximo
	//2 - Fora do range fisico da superficie - distante

	//Se está no range local de interesse - domínio físico da superfície
	if (mc(0, 0) >= surf1_knotA && mc(0, 0) < surf1_knotB && mc(1, 0) >= surf2_knotA && mc(1, 0) < surf2_knotB)
		return_value = 0;	//Houve convergência, está no range físico - forte candidato a contato
	else
	{
		if (mc(0, 0) >= (surf1_knotA - perc * surf1_range) && mc(0, 0) < (surf1_knotB + perc * surf1_range)
			&& mc(1, 0) >= (surf2_knotA - perc * surf2_range) && mc(1, 0) < (surf2_knotB + perc * surf2_range))
			return_value = 4;	//Houve convergência, está no range físico - forte candidato a contato
		else
			return_value = 2;	//Houve convergência, mas não está no range físico - deve ser monitorado com cuidado - possivelmente superfície vizinha
	}
	return return_value;
}


//Chute inicial para coordenadas convectivas do par de superfícies
void SplineElement_SplineElement::InitialGuess(SPContactData* c_data)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	SplineElement* surf1;		//Ponteiro para a superfície 1
	SplineElement* surf2;		//Ponteiro para a superfície 2
	surf1 = static_cast<SplineElement*>(db.splines[spline1_ID - 1]->sp_element[surf1_ID]);
	surf2 = static_cast<SplineElement*>(db.splines[spline2_ID - 1]->sp_element[surf2_ID]);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double tol_ortho = 1e-12;
	Matrix x_sp_AA(3);
	Matrix x_sp_BA(3);
	Matrix x_sp_AB(3);
	Matrix x_sp_BB(3);

	Matrix uAA(3);
	Matrix uBA(3);
	Matrix uAB(3);
	Matrix uBB(3);

	for (int i = 0; i < 3; i++)
	{
		uAA(i, 0) = (*surf1->d)(i, 0);
		uBA(i, 0) = (*surf1->d)(i + 6, 0);
		uAB(i, 0) = (*surf2->d)(i, 0);
		uBB(i, 0) = (*surf2->d)(i + 6, 0);
	}
	surf1->SplinePoint(surf1->knot_element[2], x_sp_AA);
	surf1->SplinePoint(surf1->knot_element[3], x_sp_BA);
	surf2->SplinePoint(surf2->knot_element[2], x_sp_AB);
	surf2->SplinePoint(surf2->knot_element[3], x_sp_BB);

	Matrix bA = 0.5*(x_sp_AA + uAA + x_sp_BA + uBA);
	Matrix tA = 0.5*(x_sp_BA + uBA - x_sp_AA - uAA);
	Matrix bB = 0.5*(x_sp_AB + uAB + x_sp_BB + uBB);
	Matrix tB = 0.5*(x_sp_BB + uBB - x_sp_AB - uAB);

	double csi_A, csi_B;

	//Estimativa dos csi's com base na configuração atual - considerando-se que são elementos retilíneos (baseado em Wriggers e Zavarise, 1997)
	csi_A = dot(bA - bB, (1.0 / (dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)))*(tB*dot(tA, tB) - tA * dot(tB, tB)));
	csi_B = -1.0*dot(bA - bB, (1.0 / (dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)))*(tA*dot(tA, tB) - tB * dot(tA, tA)));

	for (int ip = 0; ip < c_data->n_solutions; ip++)
	{
		//Preenchendo as coordenadas convectivas:
		//c_data->convective[ip][0] = 0.5*((1 - csi_A)*surf1->knot_element[2] + (1 + csi_A)*surf1->knot_element[3]);
		//c_data->convective[ip][1] = 0.5*((1 - csi_B)*surf2->knot_element[2] + (1 + csi_B)*surf2->knot_element[3]);
		c_data->convective[ip][0] = 0.5*(surf1->knot_element[2] + surf1->knot_element[3]);
		c_data->convective[ip][1] = 0.5*(surf2->knot_element[2] + surf2->knot_element[3]);
	}
}

//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
double SplineElement_SplineElement::ObjectivePhase1(Matrix& mc)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	SplineElement* surf1;		//Ponteiro para a superfície 1
	SplineElement* surf2;		//Ponteiro para a superfície 2

	double* rA;
	double* rB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	double* knotA;
	double* knotB;

	double* xAAi;
	double* xBAi;
	double* xCAi;
	double* xABi;
	double* xBBi;
	double* xCBi;

	surf1 = static_cast<SplineElement*>(db.splines[spline1_ID - 1]->sp_element[surf1_ID]);
	surf2 = static_cast<SplineElement*>(db.splines[spline2_ID - 1]->sp_element[surf2_ID]);
	rA = surf1->radius;
	rB = surf2->radius;
	dA = surf1->d->getMatrix();
	dB = surf2->d->getMatrix();
	duiA = surf1->dui->getMatrix();
	dduiA = surf1->ddui->getMatrix();
	duiB = surf2->dui->getMatrix();
	dduiB = surf2->ddui->getMatrix();
	knotA = surf1->knot_element;
	knotB = surf2->knot_element;
	xAAi = surf1->x_Ai->getMatrix();
	xBAi = surf1->x_Bi->getMatrix();
	xCAi = surf1->x_Ci->getMatrix();
	xABi = surf2->x_Ai->getMatrix();
	xBBi = surf2->x_Bi->getMatrix();
	xCBi = surf2->x_Ci->getMatrix();


	double v[356];
	double *cp = mc.getMatrix();
	double Ob;

	int b79, b85, b92, b105, b112
		, b138, b144, b151, b164
		, b171;

	if (knotA[1] != knotA[0]) {
		v[78] = 0e0;
	}
	else {
		v[78] = 0e0;
	};
	b79 = knotA[2] != knotA[1];
	if (b79) {
		v[81] = 0e0;
	}
	else {
		v[81] = 0e0;
	};
	v[82] = v[78] + v[81];
	if (b79) {
		v[84] = 0e0;
	}
	else {
		v[84] = 0e0;
	};
	b85 = knotA[3] != knotA[2];
	if (b85) {
		v[90] = 1e0 / (-knotA[2] + knotA[3]);
		v[87] = (-cp[0] + knotA[3])*v[90];
	}
	else {
		v[87] = 0e0;
	};
	v[88] = v[84] + v[87];
	if (b85) {
		v[91] = (cp[0] - knotA[2])*v[90];
	}
	else {
		v[91] = 0e0;
	};
	b92 = knotA[4] != knotA[3];
	if (b92) {
		v[94] = 0e0;
	}
	else {
		v[94] = 0e0;
	};
	v[95] = v[91] + v[94];
	if (b92) {
		v[97] = 0e0;
	}
	else {
		v[97] = 0e0;
	};
	if (knotA[5] != knotA[4]) {
		v[100] = 0e0;
	}
	else {
		v[100] = 0e0;
	};
	v[101] = v[100] + v[97];
	if (knotA[2] != knotA[0]) {
		v[104] = ((cp[0] - knotA[0])*v[82]) / (-knotA[0] + knotA[2]);
	}
	else {
		v[104] = 0e0;
	};
	b105 = knotA[3] != knotA[1];
	if (b105) {
		v[187] = v[88] / (-knotA[1] + knotA[3]);
		v[107] = (-cp[0] + knotA[3])*v[187];
	}
	else {
		v[107] = 0e0;
	};
	v[108] = v[104] + v[107];
	if (b105) {
		v[111] = (cp[0] - knotA[1])*v[187];
	}
	else {
		v[111] = 0e0;
	};
	b112 = knotA[4] != knotA[2];
	if (b112) {
		v[188] = v[95] / (-knotA[2] + knotA[4]);
		v[114] = (-cp[0] + knotA[4])*v[188];
	}
	else {
		v[114] = 0e0;
	};
	v[115] = v[111] + v[114];
	if (b112) {
		v[118] = (cp[0] - knotA[2])*v[188];
	}
	else {
		v[118] = 0e0;
	};
	if (knotA[5] != knotA[3]) {
		v[121] = ((-cp[0] + knotA[5])*v[101]) / (-knotA[3] + knotA[5]);
	}
	else {
		v[121] = 0e0;
	};
	v[122] = v[118] + v[121];
	if (knotB[1] != knotB[0]) {
		v[137] = 0e0;
	}
	else {
		v[137] = 0e0;
	};
	b138 = knotB[2] != knotB[1];
	if (b138) {
		v[140] = 0e0;
	}
	else {
		v[140] = 0e0;
	};
	v[141] = v[137] + v[140];
	if (b138) {
		v[143] = 0e0;
	}
	else {
		v[143] = 0e0;
	};
	b144 = knotB[3] != knotB[2];
	if (b144) {
		v[149] = 1e0 / (-knotB[2] + knotB[3]);
		v[146] = (-cp[1] + knotB[3])*v[149];
	}
	else {
		v[146] = 0e0;
	};
	v[147] = v[143] + v[146];
	if (b144) {
		v[150] = (cp[1] - knotB[2])*v[149];
	}
	else {
		v[150] = 0e0;
	};
	b151 = knotB[4] != knotB[3];
	if (b151) {
		v[153] = 0e0;
	}
	else {
		v[153] = 0e0;
	};
	v[154] = v[150] + v[153];
	if (b151) {
		v[156] = 0e0;
	}
	else {
		v[156] = 0e0;
	};
	if (knotB[5] != knotB[4]) {
		v[159] = 0e0;
	}
	else {
		v[159] = 0e0;
	};
	v[160] = v[156] + v[159];
	if (knotB[2] != knotB[0]) {
		v[163] = ((cp[1] - knotB[0])*v[141]) / (-knotB[0] + knotB[2]);
	}
	else {
		v[163] = 0e0;
	};
	b164 = knotB[3] != knotB[1];
	if (b164) {
		v[189] = v[147] / (-knotB[1] + knotB[3]);
		v[166] = (-cp[1] + knotB[3])*v[189];
	}
	else {
		v[166] = 0e0;
	};
	v[167] = v[163] + v[166];
	if (b164) {
		v[170] = (cp[1] - knotB[1])*v[189];
	}
	else {
		v[170] = 0e0;
	};
	b171 = knotB[4] != knotB[2];
	if (b171) {
		v[190] = v[154] / (-knotB[2] + knotB[4]);
		v[173] = (-cp[1] + knotB[4])*v[190];
	}
	else {
		v[173] = 0e0;
	};
	v[174] = v[170] + v[173];
	if (b171) {
		v[177] = (cp[1] - knotB[2])*v[190];
	}
	else {
		v[177] = 0e0;
	};
	if (knotB[5] != knotB[3]) {
		v[180] = ((-cp[1] + knotB[5])*v[160]) / (-knotB[3] + knotB[5]);
	}
	else {
		v[180] = 0e0;
	};
	v[181] = v[177] + v[180];
	(Ob) = 0.5e0*(Power(v[108] * (dA[0] + xAAi[0]) - v[167] * (dB[0] + xABi[0]) + v[115] * (dA[3] + xBAi[0]) - v[174] *
		(dB[3] + xBBi[0]) + v[122] * (dA[6] + xCAi[0]) - v[181] * (dB[6] + xCBi[0]), 2) + Power(v[108] * (dA[1] + xAAi[1])
			- v[167] * (dB[1] + xABi[1]) + v[115] * (dA[4] + xBAi[1]) - v[174] * (dB[4] + xBBi[1]) + v[122] * (dA[7] + xCAi[1])
			- v[181] * (dB[7] + xCBi[1]), 2) + Power(v[108] * (dA[2] + xAAi[2]) - v[167] * (dB[2] + xABi[2]) + v[115] * (dA[5]
				+ xBAi[2]) - v[174] * (dB[5] + xBBi[2]) + v[122] * (dA[8] + xCAi[2]) - v[181] * (dB[8] + xCBi[2]), 2));


	return Ob;
}

//Calcula o Gradiente da função objetivo - Phase 1
void SplineElement_SplineElement::GradientPhase1(Matrix& mc, Matrix& mGra)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	SplineElement* surf1;		//Ponteiro para a superfície 1
	SplineElement* surf2;		//Ponteiro para a superfície 2

	double* rA;
	double* rB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	double* knotA;
	double* knotB;

	double* xAAi;
	double* xBAi;
	double* xCAi;
	double* xABi;
	double* xBBi;
	double* xCBi;

	surf1 = static_cast<SplineElement*>(db.splines[spline1_ID - 1]->sp_element[surf1_ID]);
	surf2 = static_cast<SplineElement*>(db.splines[spline2_ID - 1]->sp_element[surf2_ID]);
	rA = surf1->radius;
	rB = surf2->radius;
	dA = surf1->d->getMatrix();
	dB = surf2->d->getMatrix();
	duiA = surf1->dui->getMatrix();
	dduiA = surf1->ddui->getMatrix();
	duiB = surf2->dui->getMatrix();
	dduiB = surf2->ddui->getMatrix();
	knotA = surf1->knot_element;
	knotB = surf2->knot_element;
	xAAi = surf1->x_Ai->getMatrix();
	xBAi = surf1->x_Bi->getMatrix();
	xCAi = surf1->x_Ci->getMatrix();
	xABi = surf2->x_Ai->getMatrix();
	xBBi = surf2->x_Bi->getMatrix();
	xCBi = surf2->x_Ci->getMatrix();

	double v[440];
	double *cp = mc.getMatrix();
	double Gra[2];

	int b79, b85, b92, b105, b112, b138, b144, b151, b164, b171;
	v[67] = dA[0] + xAAi[0];
	v[68] = dA[1] + xAAi[1];
	v[69] = dA[2] + xAAi[2];
	v[70] = dA[3] + xBAi[0];
	v[71] = dA[4] + xBAi[1];
	v[72] = dA[5] + xBAi[2];
	v[73] = dA[6] + xCAi[0];
	v[74] = dA[7] + xCAi[1];
	v[75] = dA[8] + xCAi[2];
	if (knotA[1] != knotA[0]) {
		v[78] = 0e0;
	}
	else {
		v[78] = 0e0;
	};
	b79 = knotA[2] != knotA[1];
	if (b79) {
		v[81] = 0e0;
	}
	else {
		v[81] = 0e0;
	};
	v[82] = v[78] + v[81];
	if (b79) {
		v[84] = 0e0;
	}
	else {
		v[84] = 0e0;
	};
	b85 = knotA[3] != knotA[2];
	if (b85) {
		v[90] = 1e0 / (-knotA[2] + knotA[3]);
		v[186] = -v[90];
		v[87] = (-cp[0] + knotA[3])*v[90];
	}
	else {
		v[186] = 0e0;
		v[87] = 0e0;
	};
	v[88] = v[84] + v[87];
	if (b85) {
		v[187] = v[90];
		v[91] = (cp[0] - knotA[2])*v[90];
	}
	else {
		v[187] = 0e0;
		v[91] = 0e0;
	};
	b92 = knotA[4] != knotA[3];
	if (b92) {
		v[94] = 0e0;
	}
	else {
		v[94] = 0e0;
	};
	v[95] = v[91] + v[94];
	if (b92) {
		v[97] = 0e0;
	}
	else {
		v[97] = 0e0;
	};
	if (knotA[5] != knotA[4]) {
		v[100] = 0e0;
	}
	else {
		v[100] = 0e0;
	};
	v[101] = v[100] + v[97];
	if (knotA[2] != knotA[0]) {
		v[189] = v[82] / (-knotA[0] + knotA[2]);
		v[104] = (cp[0] - knotA[0])*v[189];
	}
	else {
		v[189] = 0e0;
		v[104] = 0e0;
	};
	b105 = knotA[3] != knotA[1];
	if (b105) {
		v[232] = -cp[0] + knotA[3];
		v[110] = 1e0 / (-knotA[1] + knotA[3]);
		v[233] = v[110] * v[186];
		v[192] = -(v[110] * v[88]);
		v[190] = v[192] + v[232] * v[233];
		v[107] = -(v[192] * v[232]);
	}
	else {
		v[190] = 0e0;
		v[107] = 0e0;
	};
	v[191] = v[189] + v[190];
	v[108] = v[104] + v[107];
	if (b105) {
		v[193] = cp[0] - knotA[1];
		v[194] = -v[192] + v[193] * v[233];
		v[111] = -(v[192] * v[193]);
	}
	else {
		v[194] = 0e0;
		v[111] = 0e0;
	};
	b112 = knotA[4] != knotA[2];
	if (b112) {
		v[195] = -cp[0] + knotA[4];
		v[117] = 1e0 / (-knotA[2] + knotA[4]);
		v[234] = v[117] * v[187];
		v[198] = -(v[117] * v[95]);
		v[196] = v[198] + v[195] * v[234];
		v[114] = -(v[195] * v[198]);
	}
	else {
		v[196] = 0e0;
		v[114] = 0e0;
	};
	v[197] = v[194] + v[196];
	v[115] = v[111] + v[114];
	if (b112) {
		v[235] = cp[0] - knotA[2];
		v[199] = -v[198] + v[234] * v[235];
		v[118] = -(v[198] * v[235]);
	}
	else {
		v[199] = 0e0;
		v[118] = 0e0;
	};
	if (knotA[5] != knotA[3]) {
		v[201] = -(v[101] / (-knotA[3] + knotA[5]));
		v[121] = -((-cp[0] + knotA[5])*v[201]);
	}
	else {
		v[201] = 0e0;
		v[121] = 0e0;
	};
	v[202] = v[199] + v[201];
	v[122] = v[118] + v[121];
	v[126] = dB[0] + xABi[0];
	v[127] = dB[1] + xABi[1];
	v[128] = dB[2] + xABi[2];
	v[129] = dB[3] + xBBi[0];
	v[130] = dB[4] + xBBi[1];
	v[131] = dB[5] + xBBi[2];
	v[132] = dB[6] + xCBi[0];
	v[133] = dB[7] + xCBi[1];
	v[134] = dB[8] + xCBi[2];
	if (knotB[1] != knotB[0]) {
		v[137] = 0e0;
	}
	else {
		v[137] = 0e0;
	};
	b138 = knotB[2] != knotB[1];
	if (b138) {
		v[140] = 0e0;
	}
	else {
		v[140] = 0e0;
	};
	v[141] = v[137] + v[140];
	if (b138) {
		v[143] = 0e0;
	}
	else {
		v[143] = 0e0;
	};
	b144 = knotB[3] != knotB[2];
	if (b144) {
		v[149] = 1e0 / (-knotB[2] + knotB[3]);
		v[206] = -v[149];
		v[146] = (-cp[1] + knotB[3])*v[149];
	}
	else {
		v[206] = 0e0;
		v[146] = 0e0;
	};
	v[147] = v[143] + v[146];
	if (b144) {
		v[207] = v[149];
		v[150] = (cp[1] - knotB[2])*v[149];
	}
	else {
		v[207] = 0e0;
		v[150] = 0e0;
	};
	b151 = knotB[4] != knotB[3];
	if (b151) {
		v[153] = 0e0;
	}
	else {
		v[153] = 0e0;
	};
	v[154] = v[150] + v[153];
	if (b151) {
		v[156] = 0e0;
	}
	else {
		v[156] = 0e0;
	};
	if (knotB[5] != knotB[4]) {
		v[159] = 0e0;
	}
	else {
		v[159] = 0e0;
	};
	v[160] = v[156] + v[159];
	if (knotB[2] != knotB[0]) {
		v[209] = v[141] / (-knotB[0] + knotB[2]);
		v[163] = (cp[1] - knotB[0])*v[209];
	}
	else {
		v[209] = 0e0;
		v[163] = 0e0;
	};
	b164 = knotB[3] != knotB[1];
	if (b164) {
		v[236] = -cp[1] + knotB[3];
		v[169] = 1e0 / (-knotB[1] + knotB[3]);
		v[237] = v[169] * v[206];
		v[212] = -(v[147] * v[169]);
		v[210] = v[212] + v[236] * v[237];
		v[166] = -(v[212] * v[236]);
	}
	else {
		v[210] = 0e0;
		v[166] = 0e0;
	};
	v[211] = v[209] + v[210];
	v[167] = v[163] + v[166];
	if (b164) {
		v[213] = cp[1] - knotB[1];
		v[214] = -v[212] + v[213] * v[237];
		v[170] = -(v[212] * v[213]);
	}
	else {
		v[214] = 0e0;
		v[170] = 0e0;
	};
	b171 = knotB[4] != knotB[2];
	if (b171) {
		v[215] = -cp[1] + knotB[4];
		v[176] = 1e0 / (-knotB[2] + knotB[4]);
		v[238] = v[176] * v[207];
		v[218] = -(v[154] * v[176]);
		v[216] = v[218] + v[215] * v[238];
		v[173] = -(v[215] * v[218]);
	}
	else {
		v[216] = 0e0;
		v[173] = 0e0;
	};
	v[217] = v[214] + v[216];
	v[174] = v[170] + v[173];
	if (b171) {
		v[239] = cp[1] - knotB[2];
		v[219] = -v[218] + v[238] * v[239];
		v[177] = -(v[218] * v[239]);
	}
	else {
		v[219] = 0e0;
		v[177] = 0e0;
	};
	if (knotB[5] != knotB[3]) {
		v[221] = -(v[160] / (-knotB[3] + knotB[5]));
		v[180] = -((-cp[1] + knotB[5])*v[221]);
	}
	else {
		v[221] = 0e0;
		v[180] = 0e0;
	};
	v[222] = v[219] + v[221];
	v[181] = v[177] + v[180];
	v[240] = 2e0*(-(v[126] * v[167]) - v[129] * v[174] - v[132] * v[181] + v[108] * v[67] + v[115] * v[70] + v[122] * v[73]);
	v[241] = 2e0*(-(v[127] * v[167]) - v[130] * v[174] - v[133] * v[181] + v[108] * v[68] + v[115] * v[71] + v[122] * v[74]);
	v[242] = 2e0*(-(v[128] * v[167]) - v[131] * v[174] - v[134] * v[181] + v[108] * v[69] + v[115] * v[72] + v[122] * v[75]);
	Gra[0] = 0.5e0*(v[240] * (v[191] * v[67] + v[197] * v[70] + v[202] * v[73]) + v[241] * (v[191] * v[68] + v[197] * v[71]
		+ v[202] * v[74]) + v[242] * (v[191] * v[69] + v[197] * v[72] + v[202] * v[75]));
	Gra[1] = 0.5e0*((-(v[126] * v[211]) - v[129] * v[217] - v[132] * v[222])*v[240] + (-(v[127] * v[211]) - v[130] * v[217]
		- v[133] * v[222])*v[241] + (-(v[128] * v[211]) - v[131] * v[217] - v[134] * v[222])*v[242]);

	for (int i = 0; i < 2; i++)
		mGra(i, 0) = Gra[i];
}

//Calcula a Hessiana da função objetivo - Phase 1
void SplineElement_SplineElement::HessianPhase1(Matrix& mc, Matrix& mHes)
{
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	SplineElement* surf1;		//Ponteiro para a superfície 1
	SplineElement* surf2;		//Ponteiro para a superfície 2

	double* rA;
	double* rB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	double* knotA;
	double* knotB;

	double* xAAi;
	double* xBAi;
	double* xCAi;
	double* xABi;
	double* xBBi;
	double* xCBi;

	surf1 = static_cast<SplineElement*>(db.splines[spline1_ID - 1]->sp_element[surf1_ID]);
	surf2 = static_cast<SplineElement*>(db.splines[spline2_ID - 1]->sp_element[surf2_ID]);

	rA = surf1->radius;
	rB = surf2->radius;
	dA = surf1->d->getMatrix();
	dB = surf2->d->getMatrix();
	duiA = surf1->dui->getMatrix();
	dduiA = surf1->ddui->getMatrix();
	duiB = surf2->dui->getMatrix();
	dduiB = surf2->ddui->getMatrix();
	knotA = surf1->knot_element;
	knotB = surf2->knot_element;
	xAAi = surf1->x_Ai->getMatrix();
	xBAi = surf1->x_Bi->getMatrix();
	xCAi = surf1->x_Ci->getMatrix();
	xABi = surf2->x_Ai->getMatrix();
	xBBi = surf2->x_Bi->getMatrix();
	xCBi = surf2->x_Ci->getMatrix();

	double v[554];
	double *cp = mc.getMatrix();
	double Hes[2][2];

	int  b79, b85, b92, b105, b112
		, b138, b144, b151, b164
		, b171;
	v[67] = dA[0] + xAAi[0];
	v[68] = dA[1] + xAAi[1];
	v[69] = dA[2] + xAAi[2];
	v[70] = dA[3] + xBAi[0];
	v[71] = dA[4] + xBAi[1];
	v[72] = dA[5] + xBAi[2];
	v[73] = dA[6] + xCAi[0];
	v[74] = dA[7] + xCAi[1];
	v[75] = dA[8] + xCAi[2];
	if (knotA[1] != knotA[0]) {
		v[78] = 0e0;
	}
	else {
		v[78] = 0e0;
	};
	b79 = knotA[2] != knotA[1];
	if (b79) {
		v[81] = 0e0;
	}
	else {
		v[81] = 0e0;
	};
	v[82] = v[78] + v[81];
	if (b79) {
		v[84] = 0e0;
	}
	else {
		v[84] = 0e0;
	};
	b85 = knotA[3] != knotA[2];
	if (b85) {
		v[90] = 1e0 / (-knotA[2] + knotA[3]);
		v[186] = -v[90];
		v[231] = -v[90];
		v[87] = (-cp[0] + knotA[3])*v[90];
	}
	else {
		v[186] = 0e0;
		v[231] = 0e0;
		v[87] = 0e0;
	};
	v[88] = v[84] + v[87];
	if (b85) {
		v[187] = v[90];
		v[232] = v[90];
		v[91] = (cp[0] - knotA[2])*v[90];
	}
	else {
		v[187] = 0e0;
		v[232] = 0e0;
		v[91] = 0e0;
	};
	b92 = knotA[4] != knotA[3];
	if (b92) {
		v[94] = 0e0;
	}
	else {
		v[94] = 0e0;
	};
	v[95] = v[91] + v[94];
	if (b92) {
		v[97] = 0e0;
	}
	else {
		v[97] = 0e0;
	};
	if (knotA[5] != knotA[4]) {
		v[100] = 0e0;
	}
	else {
		v[100] = 0e0;
	};
	v[101] = v[100] + v[97];
	if (knotA[2] != knotA[0]) {
		v[233] = v[82] / (-knotA[0] + knotA[2]);
		v[189] = v[233];
		v[234] = v[233];
		v[104] = (cp[0] - knotA[0])*v[233];
	}
	else {
		v[189] = 0e0;
		v[234] = 0e0;
		v[104] = 0e0;
	};
	b105 = knotA[3] != knotA[1];
	if (b105) {
		v[237] = -cp[0] + knotA[3];
		v[110] = 1e0 / (-knotA[1] + knotA[3]);
		v[296] = v[110] * v[186];
		v[235] = -(v[110] * v[231]);
		v[240] = v[235] - v[296];
		v[192] = -(v[110] * v[88]);
		v[236] = v[240];
		v[190] = v[192] + v[237] * v[296];
		v[238] = v[192] - v[235] * v[237];
		v[107] = -(v[192] * v[237]);
	}
	else {
		v[236] = 0e0;
		v[190] = 0e0;
		v[238] = 0e0;
		v[107] = 0e0;
	};
	v[239] = v[234] + v[238];
	v[191] = v[189] + v[190];
	v[108] = v[104] + v[107];
	if (b105) {
		v[193] = cp[0] - knotA[1];
		v[241] = -v[240];
		v[194] = -v[192] + v[193] * v[296];
		v[242] = -v[192] - v[193] * v[235];
		v[111] = -(v[192] * v[193]);
	}
	else {
		v[241] = 0e0;
		v[194] = 0e0;
		v[242] = 0e0;
		v[111] = 0e0;
	};
	b112 = knotA[4] != knotA[2];
	if (b112) {
		v[195] = -cp[0] + knotA[4];
		v[117] = 1e0 / (-knotA[2] + knotA[4]);
		v[297] = v[117] * v[187];
		v[243] = -(v[117] * v[232]);
		v[248] = v[243] - v[297];
		v[198] = -(v[117] * v[95]);
		v[244] = v[248];
		v[196] = v[198] + v[195] * v[297];
		v[245] = v[198] - v[195] * v[243];
		v[114] = -(v[195] * v[198]);
	}
	else {
		v[244] = 0e0;
		v[196] = 0e0;
		v[245] = 0e0;
		v[114] = 0e0;
	};
	v[247] = v[242] + v[245];
	v[246] = v[241] + v[244];
	v[197] = v[194] + v[196];
	v[115] = v[111] + v[114];
	if (b112) {
		v[250] = cp[0] - knotA[2];
		v[249] = -v[248];
		v[199] = -v[198] + v[250] * v[297];
		v[251] = -v[198] - v[243] * v[250];
		v[118] = -(v[198] * v[250]);
	}
	else {
		v[249] = 0e0;
		v[199] = 0e0;
		v[251] = 0e0;
		v[118] = 0e0;
	};
	if (knotA[5] != knotA[3]) {
		v[252] = -(v[101] / (-knotA[3] + knotA[5]));
		v[201] = v[252];
		v[253] = v[252];
		v[121] = -((-cp[0] + knotA[5])*v[252]);
	}
	else {
		v[201] = 0e0;
		v[253] = 0e0;
		v[121] = 0e0;
	};
	v[257] = v[251] + v[253];
	v[305] = 2e0*(v[239] * v[69] + v[247] * v[72] + v[257] * v[75]);
	v[304] = 2e0*(v[239] * v[68] + v[247] * v[71] + v[257] * v[74]);
	v[303] = 2e0*(v[239] * v[67] + v[247] * v[70] + v[257] * v[73]);
	v[202] = v[199] + v[201];
	v[205] = v[191] * v[69] + v[197] * v[72] + v[202] * v[75];
	v[204] = v[191] * v[68] + v[197] * v[71] + v[202] * v[74];
	v[203] = v[191] * v[67] + v[197] * v[70] + v[202] * v[73];
	v[122] = v[118] + v[121];
	v[126] = dB[0] + xABi[0];
	v[127] = dB[1] + xABi[1];
	v[128] = dB[2] + xABi[2];
	v[129] = dB[3] + xBBi[0];
	v[130] = dB[4] + xBBi[1];
	v[131] = dB[5] + xBBi[2];
	v[132] = dB[6] + xCBi[0];
	v[133] = dB[7] + xCBi[1];
	v[134] = dB[8] + xCBi[2];
	if (knotB[1] != knotB[0]) {
		v[137] = 0e0;
	}
	else {
		v[137] = 0e0;
	};
	b138 = knotB[2] != knotB[1];
	if (b138) {
		v[140] = 0e0;
	}
	else {
		v[140] = 0e0;
	};
	v[141] = v[137] + v[140];
	if (b138) {
		v[143] = 0e0;
	}
	else {
		v[143] = 0e0;
	};
	b144 = knotB[3] != knotB[2];
	if (b144) {
		v[149] = 1e0 / (-knotB[2] + knotB[3]);
		v[206] = -v[149];
		v[261] = -v[149];
		v[146] = (-cp[1] + knotB[3])*v[149];
	}
	else {
		v[206] = 0e0;
		v[261] = 0e0;
		v[146] = 0e0;
	};
	v[147] = v[143] + v[146];
	if (b144) {
		v[207] = v[149];
		v[262] = v[149];
		v[150] = (cp[1] - knotB[2])*v[149];
	}
	else {
		v[207] = 0e0;
		v[262] = 0e0;
		v[150] = 0e0;
	};
	b151 = knotB[4] != knotB[3];
	if (b151) {
		v[153] = 0e0;
	}
	else {
		v[153] = 0e0;
	};
	v[154] = v[150] + v[153];
	if (b151) {
		v[156] = 0e0;
	}
	else {
		v[156] = 0e0;
	};
	if (knotB[5] != knotB[4]) {
		v[159] = 0e0;
	}
	else {
		v[159] = 0e0;
	};
	v[160] = v[156] + v[159];
	if (knotB[2] != knotB[0]) {
		v[263] = v[141] / (-knotB[0] + knotB[2]);
		v[209] = v[263];
		v[264] = v[263];
		v[163] = (cp[1] - knotB[0])*v[263];
	}
	else {
		v[209] = 0e0;
		v[264] = 0e0;
		v[163] = 0e0;
	};
	b164 = knotB[3] != knotB[1];
	if (b164) {
		v[267] = -cp[1] + knotB[3];
		v[169] = 1e0 / (-knotB[1] + knotB[3]);
		v[298] = v[169] * v[206];
		v[265] = -(v[169] * v[261]);
		v[270] = v[265] - v[298];
		v[212] = -(v[147] * v[169]);
		v[266] = v[270];
		v[210] = v[212] + v[267] * v[298];
		v[268] = v[212] - v[265] * v[267];
		v[166] = -(v[212] * v[267]);
	}
	else {
		v[266] = 0e0;
		v[210] = 0e0;
		v[268] = 0e0;
		v[166] = 0e0;
	};
	v[269] = v[264] + v[268];
	v[211] = v[209] + v[210];
	v[167] = v[163] + v[166];
	if (b164) {
		v[213] = cp[1] - knotB[1];
		v[271] = -v[270];
		v[214] = -v[212] + v[213] * v[298];
		v[272] = -v[212] - v[213] * v[265];
		v[170] = -(v[212] * v[213]);
	}
	else {
		v[271] = 0e0;
		v[214] = 0e0;
		v[272] = 0e0;
		v[170] = 0e0;
	};
	b171 = knotB[4] != knotB[2];
	if (b171) {
		v[215] = -cp[1] + knotB[4];
		v[176] = 1e0 / (-knotB[2] + knotB[4]);
		v[299] = v[176] * v[207];
		v[273] = -(v[176] * v[262]);
		v[278] = v[273] - v[299];
		v[218] = -(v[154] * v[176]);
		v[274] = v[278];
		v[216] = v[218] + v[215] * v[299];
		v[275] = v[218] - v[215] * v[273];
		v[173] = -(v[215] * v[218]);
	}
	else {
		v[274] = 0e0;
		v[216] = 0e0;
		v[275] = 0e0;
		v[173] = 0e0;
	};
	v[277] = v[272] + v[275];
	v[276] = v[271] + v[274];
	v[217] = v[214] + v[216];
	v[174] = v[170] + v[173];
	if (b171) {
		v[280] = cp[1] - knotB[2];
		v[279] = -v[278];
		v[219] = -v[218] + v[280] * v[299];
		v[281] = -v[218] - v[273] * v[280];
		v[177] = -(v[218] * v[280]);
	}
	else {
		v[279] = 0e0;
		v[219] = 0e0;
		v[281] = 0e0;
		v[177] = 0e0;
	};
	if (knotB[5] != knotB[3]) {
		v[282] = -(v[160] / (-knotB[3] + knotB[5]));
		v[221] = v[282];
		v[283] = v[282];
		v[180] = -((-cp[1] + knotB[5])*v[282]);
	}
	else {
		v[221] = 0e0;
		v[283] = 0e0;
		v[180] = 0e0;
	};
	v[287] = v[281] + v[283];
	v[308] = 2e0*(v[128] * v[269] + v[131] * v[277] + v[134] * v[287]);
	v[307] = 2e0*(v[127] * v[269] + v[130] * v[277] + v[133] * v[287]);
	v[306] = 2e0*(v[126] * v[269] + v[129] * v[277] + v[132] * v[287]);
	v[222] = v[219] + v[221];
	v[225] = v[128] * v[211] + v[131] * v[217] + v[134] * v[222];
	v[224] = v[127] * v[211] + v[130] * v[217] + v[133] * v[222];
	v[223] = v[126] * v[211] + v[129] * v[217] + v[132] * v[222];
	v[181] = v[177] + v[180];
	v[300] = 2e0*(-(v[126] * v[167]) - v[129] * v[174] - v[132] * v[181] + v[108] * v[67] + v[115] * v[70] + v[122] * v[73]);
	v[301] = 2e0*(-(v[127] * v[167]) - v[130] * v[174] - v[133] * v[181] + v[108] * v[68] + v[115] * v[71] + v[122] * v[74]);
	v[302] = 2e0*(-(v[128] * v[167]) - v[131] * v[174] - v[134] * v[181] + v[108] * v[69] + v[115] * v[72] + v[122] * v[75]);
	Hes[0][0] = 0.5e0*(v[203] * v[303] + v[204] * v[304] + v[205] * v[305] + v[300] * (v[236] * v[67] + v[246] * v[70]
		+ v[249] * v[73]) + v[301] * (v[236] * v[68] + v[246] * v[71] + v[249] * v[74]) + v[302] * (v[236] * v[69] + v[246] * v[72]
			+ v[249] * v[75]));
	Hes[0][1] = 0.5e0*(-(v[203] * v[306]) - v[204] * v[307] - v[205] * v[308]);
	Hes[1][0] = 0.5e0*(-(v[223] * v[303]) - v[224] * v[304] - v[225] * v[305]);
	Hes[1][1] = 0.5e0*((-(v[126] * v[266]) - v[129] * v[276] - v[132] * v[279])*v[300] + (-(v[127] * v[266])
		- v[130] * v[276] - v[133] * v[279])*v[301] + (-(v[128] * v[266]) - v[131] * v[276] - v[134] * v[279])*v[302]
		+ v[223] * v[306] + v[224] * v[307] + v[225] * v[308]);

	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			mHes(i, j) = Hes[i][j];

}

void SplineElement_SplineElement::ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn_n, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft)
{
	double v[4607];		//variável temporária - AceGen
	//Zerando matrizes e vetores
	for (int i = 0; i < 18; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < 18; j++)
			Kc[i][j] = 0.0;
	}
	///////////////////////////////Ponteiros e variáveis para facilitar acesso//////////////////////////////////////////
	SplineElement* surf1;		//Ponteiro para a superfície 1
	SplineElement* surf2;		//Ponteiro para a superfície 2
	double* rA;
	double* rB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	double* knotA;
	double* knotB;
	double* xAAi;
	double* xBAi;
	double* xCAi;
	double* xABi;
	double* xBBi;
	double* xCBi;
	surf1 = static_cast<SplineElement*>(db.splines[spline1_ID - 1]->sp_element[surf1_ID]);
	surf2 = static_cast<SplineElement*>(db.splines[spline2_ID - 1]->sp_element[surf2_ID]);

	rA = surf1->radius;
	rB = surf2->radius;

	dA = surf1->d->getMatrix();
	dB = surf2->d->getMatrix();

	duiA = surf1->dui->getMatrix();
	dduiA = surf1->ddui->getMatrix();

	duiB = surf2->dui->getMatrix();
	dduiB = surf2->ddui->getMatrix();

	knotA = surf1->knot_element;
	knotB = surf2->knot_element;

	xAAi = surf1->x_Ai->getMatrix();
	xBAi = surf1->x_Bi->getMatrix();
	xCAi = surf1->x_Ci->getMatrix();

	xABi = surf2->x_Ai->getMatrix();
	xBBi = surf2->x_Bi->getMatrix();
	xCBi = surf2->x_Ci->getMatrix();

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

	//double epsn_n_v = 2.4274;
	//double *epsn_n = &epsn_n_v;

#pragma region AceGen
	double v01; double v010; double v011; double v012; double v013; double v014;
	double v015; double v016; double v017; double v018; double v019; double v02;
	double v020; double v021; double v022; double v023; double v024; double v025;
	double v026; double v027; double v03; double v04; double v05; double v06; double v07;
	double v08; double v09;
	int i667, i792, i1213, i1513, i2194, i2267, i2507, i2508, i2509, i2510, i2511, i2512, i2513
		, i2514, i2515, i2516, i2517, i2518, i2519, i2520, i2521, i2522, i2523, i2524, i2587, i2588
		, i2589, i2590, i2591, i2592, i2593, i2594, i2595, i2596, i2597, i2598, i2599, i2600, i2601
		, i2602, i2603, i2604, i2652, i2653, i2654, i2655, i2656, i2657, i2658, i2659, i2660, i2661
		, i2662, i2663, i2664, i2665, i2666, i2667, i2668, i2669, b22, b23, b128, b131
		, b137, b144, b150, b154, b157, b164
		, b171
		, b230, b233, b239, b246, b252
		, b256, b259, b266, b273
		, b480, b598
		, b632, b722, b724, b726, b728, b729, b731, b733, b735, b737, b739
		, b741, b743, b744, b746, b748, b750, b815, b818, b822, b825, b828, b831, b834, b835, b836
		, b839, b843, b846, b849, b852, b855, b856, b1089, b1091, b1092, b1094, b1095, b1096, b1098
		, b1100, b1102, b1104, b1105, b1107, b1108, b1109, b1111, b1113, b1218, b1219, b1220
		, b1358, b1389, b1460, b1463, b1464, b1467, b1469, b1470, b1472, b1474, b1486, b1489, b1490
		, b1493, b1495, b1496, b1498, b1500, b1618, b1664, b1776, b1777, b1793, b2027, b2046, b2099
		, b2102, b2103, b2106, b2109, b2110, b2112, b2114, b2129, b2132, b2133, b2136, b2139, b2140
		, b2142, b2144, b2203, b2205, b2206, b2208, b2209, b2210, b2212, b2214, b2216, b2218, b2219
		, b2221, b2222, b2223, b2225, b2227, b2289, b2292, b2296, b2299, b2302, b2305, b2307, b2308
		, b2309, b2312, b2316, b2319, b2322, b2325, b2327, b2328, b2389, b2391, b2392, b2394, b2395
		, b2396, b2398, b2400, b2401, b2403, b2404, b2406, b2407, b2408, b2410, b2412;
	v[1] = cp[0];
	v[2] = cp[1];
	v[3] = ci[0];
	v[4] = ci[1];
	v[5] = gti[0];
	v[6] = gti[1];
	v[7] = gti[2];
	v[9] = (*epsn_n);
	v[10] = (*epst);
	v[11] = (*mus);
	v[12] = (*mud);
	v[13] = (*cn);
	v[14] = (*ct);
	v[15] = (*a4);
	v[2586] = v[14] * v[15];
	v[16] = (*a5);
	v[17] = (*a6);
	v[18] = invH[0][0];
	v[19] = invH[0][1];
	v[20] = invH[1][0];
	v[21] = invH[1][1];
	b22 = (*stick);
	b23 = (*previouscontact);
	v[24] = dA[0];
	v[25] = dA[1];
	v[26] = dA[2];
	v[27] = dA[3];
	v[28] = dA[4];
	v[29] = dA[5];
	v[30] = dA[6];
	v[31] = dA[7];
	v[32] = dA[8];
	v[52] = xAAi[0];
	v[53] = xAAi[1];
	v[54] = xAAi[2];
	v[55] = xBAi[0];
	v[56] = xBAi[1];
	v[57] = xBAi[2];
	v[58] = xCAi[0];
	v[59] = xCAi[1];
	v[60] = xCAi[2];
	v[61] = knotA[0];
	v[62] = knotA[1];
	v[63] = knotA[2];
	v[64] = knotA[3];
	v[65] = knotA[4];
	v[66] = knotA[5];
	v[67] = dB[0];
	v[68] = dB[1];
	v[69] = dB[2];
	v[70] = dB[3];
	v[71] = dB[4];
	v[72] = dB[5];
	v[73] = dB[6];
	v[74] = dB[7];
	v[75] = dB[8];
	v[95] = xABi[0];
	v[96] = xABi[1];
	v[97] = xABi[2];
	v[98] = xBBi[0];
	v[99] = xBBi[1];
	v[100] = xBBi[2];
	v[101] = xCBi[0];
	v[102] = xCBi[1];
	v[103] = xCBi[2];
	v[104] = knotB[0];
	v[105] = knotB[1];
	v[106] = knotB[2];
	v[107] = knotB[3];
	v[108] = knotB[4];
	v[109] = knotB[5];
	v[119] = v[24] + v[52];
	v[120] = v[25] + v[53];
	v[121] = v[26] + v[54];
	v[122] = v[27] + v[55];
	v[123] = v[28] + v[56];
	v[124] = v[29] + v[57];
	v[125] = v[30] + v[58];
	v[126] = v[31] + v[59];
	v[127] = v[32] + v[60];
	b128 = v[62] != v[61];
	if (b128) {
		v[130] = 0e0;
	}
	else {
		v[130] = 0e0;
	};
	b131 = v[63] != v[62];
	if (b131) {
		v[133] = 0e0;
	}
	else {
		v[133] = 0e0;
	};
	v[134] = v[130] + v[133];
	if (b131) {
		v[136] = 0e0;
	}
	else {
		v[136] = 0e0;
	};
	b137 = v[64] != v[63];
	if (b137) {
		v[142] = 1e0 / (-v[63] + v[64]);
		v[139] = v[142] * (-v[3] + v[64]);
	}
	else {
		v[139] = 0e0;
	};
	v[140] = v[136] + v[139];
	if (b137) {
		v[143] = v[142] * (v[3] - v[63]);
	}
	else {
		v[143] = 0e0;
	};
	b144 = v[65] != v[64];
	if (b144) {
		v[146] = 0e0;
	}
	else {
		v[146] = 0e0;
	};
	v[147] = v[143] + v[146];
	if (b144) {
		v[149] = 0e0;
	}
	else {
		v[149] = 0e0;
	};
	b150 = v[66] != v[65];
	if (b150) {
		v[152] = 0e0;
	}
	else {
		v[152] = 0e0;
	};
	v[153] = v[149] + v[152];
	b154 = v[63] != v[61];
	if (b154) {
		v[196] = 1e0 / (-v[61] + v[63]);
		v[156] = v[134] * v[196] * (v[3] - v[61]);
	}
	else {
		v[156] = 0e0;
	};
	b157 = v[64] != v[62];
	if (b157) {
		v[162] = 1e0 / (-v[62] + v[64]);
		v[2454] = v[140] * v[162];
		v[159] = v[2454] * (-v[3] + v[64]);
	}
	else {
		v[159] = 0e0;
	};
	v[160] = v[156] + v[159];
	if (b157) {
		v[163] = v[2454] * (v[3] - v[62]);
	}
	else {
		v[163] = 0e0;
	};
	b164 = v[65] != v[63];
	if (b164) {
		v[169] = 1e0 / (-v[63] + v[65]);
		v[2455] = v[147] * v[169];
		v[166] = v[2455] * (-v[3] + v[65]);
	}
	else {
		v[166] = 0e0;
	};
	v[167] = v[163] + v[166];
	if (b164) {
		v[170] = v[2455] * (v[3] - v[63]);
	}
	else {
		v[170] = 0e0;
	};
	b171 = v[66] != v[64];
	if (b171) {
		v[209] = 1e0 / (-v[64] + v[66]);
		v[173] = v[153] * v[209] * (-v[3] + v[66]);
	}
	else {
		v[173] = 0e0;
	};
	v[174] = v[170] + v[173];
	if (b128) {
		v[176] = 0e0;
	}
	else {
		v[176] = 0e0;
	};
	if (b131) {
		v[178] = 0e0;
	}
	else {
		v[178] = 0e0;
	};
	v[179] = v[176] + v[178];
	if (b131) {
		v[181] = 0e0;
	}
	else {
		v[181] = 0e0;
	};
	if (b137) {
		v[370] = -v[142];
		v[538] = -v[142];
		v[1162] = -v[142];
		v[183] = v[142] * (-v[1] + v[64]);
	}
	else {
		v[370] = 0e0;
		v[538] = 0e0;
		v[1162] = 0e0;
		v[183] = 0e0;
	};
	v[184] = v[181] + v[183];
	if (b137) {
		v[371] = v[142];
		v[539] = v[142];
		v[1163] = v[142];
		v[186] = v[142] * (v[1] - v[63]);
	}
	else {
		v[371] = 0e0;
		v[539] = 0e0;
		v[1163] = 0e0;
		v[186] = 0e0;
	};
	if (b144) {
		v[188] = 0e0;
	}
	else {
		v[188] = 0e0;
	};
	v[189] = v[186] + v[188];
	if (b144) {
		v[191] = 0e0;
	}
	else {
		v[191] = 0e0;
	};
	if (b150) {
		v[193] = 0e0;
	}
	else {
		v[193] = 0e0;
	};
	v[194] = v[191] + v[193];
	if (b154) {
		v[540] = v[179] * v[196];
		v[372] = v[540];
		v[541] = v[540];
		v[1165] = v[540];
		v[197] = v[540] * (v[1] - v[61]);
	}
	else {
		v[372] = 0e0;
		v[541] = 0e0;
		v[1165] = 0e0;
		v[197] = 0e0;
	};
	if (b157) {
		v[745] = -v[1] + v[64];
		v[2456] = v[162] * v[745];
		v[375] = -(v[162] * v[184]);
		v[373] = v[2456] * v[370] + v[375];
		v[543] = v[375] + v[2456] * v[538];
		v[1166] = v[1162] * v[2456] + v[375];
		v[199] = v[184] * v[2456];
	}
	else {
		v[373] = 0e0;
		v[543] = 0e0;
		v[1166] = 0e0;
		v[199] = 0e0;
	};
	v[1167] = v[1165] + v[1166];
	v[544] = v[541] + v[543];
	v[374] = v[372] + v[373];
	v[200] = v[197] + v[199];
	v[2506] = v[15] * v[200];
	if (b157) {
		v[2457] = v[162] * (v[1] - v[62]);
		v[377] = v[2457] * v[370] - v[375];
		v[545] = -v[375] + v[2457] * v[538];
		v[1168] = v[1162] * v[2457] - v[375];
		v[202] = v[184] * v[2457];
	}
	else {
		v[377] = 0e0;
		v[545] = 0e0;
		v[1168] = 0e0;
		v[202] = 0e0;
	};
	if (b164) {
		v[381] = -(v[169] * v[189]);
		v[378] = -v[1] + v[65];
		v[2458] = v[169] * v[378];
		v[379] = v[2458] * v[371] + v[381];
		v[546] = v[381] + v[2458] * v[539];
		v[1169] = v[1163] * v[2458] + v[381];
		v[204] = v[189] * v[2458];
	}
	else {
		v[379] = 0e0;
		v[546] = 0e0;
		v[1169] = 0e0;
		v[204] = 0e0;
	};
	v[1170] = v[1168] + v[1169];
	v[547] = v[545] + v[546];
	v[380] = v[377] + v[379];
	v[205] = v[202] + v[204];
	v[2505] = v[15] * v[205];
	if (b164) {
		v[2459] = v[169] * (v[1] - v[63]);
		v[382] = v[2459] * v[371] - v[381];
		v[549] = -v[381] + v[2459] * v[539];
		v[1171] = v[1163] * v[2459] - v[381];
		v[207] = v[189] * v[2459];
	}
	else {
		v[382] = 0e0;
		v[549] = 0e0;
		v[1171] = 0e0;
		v[207] = 0e0;
	};
	if (b171) {
		v[550] = -(v[194] * v[209]);
		v[383] = v[550];
		v[551] = v[550];
		v[1173] = v[550];
		v[210] = -(v[550] * (-v[1] + v[66]));
	}
	else {
		v[383] = 0e0;
		v[551] = 0e0;
		v[1173] = 0e0;
		v[210] = 0e0;
	};
	v[1174] = v[1171] + v[1173];
	v[1177] = v[1167] * v[121] + v[1170] * v[124] + v[1174] * v[127];
	v[1176] = v[1167] * v[120] + v[1170] * v[123] + v[1174] * v[126];
	v[1175] = v[1167] * v[119] + v[1170] * v[122] + v[1174] * v[125];
	v[552] = v[549] + v[551];
	v[555] = v[121] * v[544] + v[124] * v[547] + v[127] * v[552];
	v[554] = v[120] * v[544] + v[123] * v[547] + v[126] * v[552];
	v[553] = v[119] * v[544] + v[122] * v[547] + v[125] * v[552];
	v[384] = v[382] + v[383];
	v[387] = v[121] * v[374] + v[124] * v[380] + v[127] * v[384];
	v[386] = v[120] * v[374] + v[123] * v[380] + v[126] * v[384];
	v[385] = v[119] * v[374] + v[122] * v[380] + v[125] * v[384];
	v[211] = v[207] + v[210];
	v[2504] = v[15] * v[211];
	v[221] = v[67] + v[95];
	v[222] = v[68] + v[96];
	v[223] = v[69] + v[97];
	v[224] = v[70] + v[98];
	v[225] = v[71] + v[99];
	v[226] = v[100] + v[72];
	v[227] = v[101] + v[73];
	v[228] = v[102] + v[74];
	v[229] = v[103] + v[75];
	b230 = v[105] != v[104];
	if (b230) {
		v[232] = 0e0;
	}
	else {
		v[232] = 0e0;
	};
	b233 = v[106] != v[105];
	if (b233) {
		v[235] = 0e0;
	}
	else {
		v[235] = 0e0;
	};
	v[236] = v[232] + v[235];
	if (b233) {
		v[238] = 0e0;
	}
	else {
		v[238] = 0e0;
	};
	b239 = v[107] != v[106];
	if (b239) {
		v[244] = 1e0 / (-v[106] + v[107]);
		v[241] = v[244] * (v[107] - v[4]);
	}
	else {
		v[241] = 0e0;
	};
	v[242] = v[238] + v[241];
	if (b239) {
		v[245] = v[244] * (-v[106] + v[4]);
	}
	else {
		v[245] = 0e0;
	};
	b246 = v[108] != v[107];
	if (b246) {
		v[248] = 0e0;
	}
	else {
		v[248] = 0e0;
	};
	v[249] = v[245] + v[248];
	if (b246) {
		v[251] = 0e0;
	}
	else {
		v[251] = 0e0;
	};
	b252 = v[109] != v[108];
	if (b252) {
		v[254] = 0e0;
	}
	else {
		v[254] = 0e0;
	};
	v[255] = v[251] + v[254];
	b256 = v[106] != v[104];
	if (b256) {
		v[298] = 1e0 / (-v[104] + v[106]);
		v[258] = v[236] * v[298] * (-v[104] + v[4]);
	}
	else {
		v[258] = 0e0;
	};
	b259 = v[107] != v[105];
	if (b259) {
		v[264] = 1e0 / (-v[105] + v[107]);
		v[2460] = v[242] * v[264];
		v[261] = v[2460] * (v[107] - v[4]);
	}
	else {
		v[261] = 0e0;
	};
	v[262] = v[258] + v[261];
	if (b259) {
		v[265] = v[2460] * (-v[105] + v[4]);
	}
	else {
		v[265] = 0e0;
	};
	b266 = v[108] != v[106];
	if (b266) {
		v[271] = 1e0 / (-v[106] + v[108]);
		v[2461] = v[249] * v[271];
		v[268] = v[2461] * (v[108] - v[4]);
	}
	else {
		v[268] = 0e0;
	};
	v[269] = v[265] + v[268];
	if (b266) {
		v[272] = v[2461] * (-v[106] + v[4]);
	}
	else {
		v[272] = 0e0;
	};
	b273 = v[109] != v[107];
	if (b273) {
		v[311] = 1e0 / (-v[107] + v[109]);
		v[275] = v[255] * v[311] * (v[109] - v[4]);
	}
	else {
		v[275] = 0e0;
	};
	v[276] = v[272] + v[275];
	if (b230) {
		v[278] = 0e0;
	}
	else {
		v[278] = 0e0;
	};
	if (b233) {
		v[280] = 0e0;
	}
	else {
		v[280] = 0e0;
	};
	v[281] = v[278] + v[280];
	if (b233) {
		v[283] = 0e0;
	}
	else {
		v[283] = 0e0;
	};
	if (b239) {
		v[389] = -v[244];
		v[556] = -v[244];
		v[1178] = -v[244];
		v[285] = (v[107] - v[2]) * v[244];
	}
	else {
		v[389] = 0e0;
		v[556] = 0e0;
		v[1178] = 0e0;
		v[285] = 0e0;
	};
	v[286] = v[283] + v[285];
	if (b239) {
		v[390] = v[244];
		v[557] = v[244];
		v[1179] = v[244];
		v[288] = (-v[106] + v[2]) * v[244];
	}
	else {
		v[390] = 0e0;
		v[557] = 0e0;
		v[1179] = 0e0;
		v[288] = 0e0;
	};
	if (b246) {
		v[290] = 0e0;
	}
	else {
		v[290] = 0e0;
	};
	v[291] = v[288] + v[290];
	if (b246) {
		v[293] = 0e0;
	}
	else {
		v[293] = 0e0;
	};
	if (b252) {
		v[295] = 0e0;
	}
	else {
		v[295] = 0e0;
	};
	v[296] = v[293] + v[295];
	if (b256) {
		v[558] = v[281] * v[298];
		v[391] = v[558];
		v[559] = v[558];
		v[1181] = v[558];
		v[299] = (-v[104] + v[2]) * v[558];
	}
	else {
		v[391] = 0e0;
		v[559] = 0e0;
		v[1181] = 0e0;
		v[299] = 0e0;
	};
	if (b259) {
		v[730] = v[107] - v[2];
		v[2462] = v[264] * v[730];
		v[394] = -(v[264] * v[286]);
		v[392] = v[2462] * v[389] + v[394];
		v[561] = v[394] + v[2462] * v[556];
		v[1182] = v[1178] * v[2462] + v[394];
		v[301] = v[2462] * v[286];
	}
	else {
		v[392] = 0e0;
		v[561] = 0e0;
		v[1182] = 0e0;
		v[301] = 0e0;
	};
	v[1183] = v[1181] + v[1182];
	v[562] = v[559] + v[561];
	v[393] = v[391] + v[392];
	v[302] = v[299] + v[301];
	v[2503] = v[15] * v[302];
	if (b259) {
		v[2463] = (-v[105] + v[2]) * v[264];
		v[396] = v[2463] * v[389] - v[394];
		v[563] = -v[394] + v[2463] * v[556];
		v[1184] = v[1178] * v[2463] - v[394];
		v[304] = v[2463] * v[286];
	}
	else {
		v[396] = 0e0;
		v[563] = 0e0;
		v[1184] = 0e0;
		v[304] = 0e0;
	};
	if (b266) {
		v[400] = -(v[271] * v[291]);
		v[397] = v[108] - v[2];
		v[2464] = v[271] * v[397];
		v[398] = v[2464] * v[390] + v[400];
		v[564] = v[400] + v[2464] * v[557];
		v[1185] = v[1179] * v[2464] + v[400];
		v[306] = v[2464] * v[291];
	}
	else {
		v[398] = 0e0;
		v[564] = 0e0;
		v[1185] = 0e0;
		v[306] = 0e0;
	};
	v[1186] = v[1184] + v[1185];
	v[565] = v[563] + v[564];
	v[399] = v[396] + v[398];
	v[307] = v[304] + v[306];
	v[2502] = v[15] * v[307];
	if (b266) {
		v[2465] = (-v[106] + v[2]) * v[271];
		v[401] = v[2465] * v[390] - v[400];
		v[567] = -v[400] + v[2465] * v[557];
		v[1187] = v[1179] * v[2465] - v[400];
		v[309] = v[2465] * v[291];
	}
	else {
		v[401] = 0e0;
		v[567] = 0e0;
		v[1187] = 0e0;
		v[309] = 0e0;
	};
	if (b273) {
		v[568] = -(v[296] * v[311]);
		v[402] = v[568];
		v[569] = v[568];
		v[1189] = v[568];
		v[312] = -((v[109] - v[2]) * v[568]);
	}
	else {
		v[402] = 0e0;
		v[569] = 0e0;
		v[1189] = 0e0;
		v[312] = 0e0;
	};
	v[1190] = v[1187] + v[1189];
	v[1193] = v[1183] * v[223] + v[1186] * v[226] + v[1190] * v[229];
	v[1192] = v[1183] * v[222] + v[1186] * v[225] + v[1190] * v[228];
	v[1191] = v[1183] * v[221] + v[1186] * v[224] + v[1190] * v[227];
	v[570] = v[567] + v[569];
	v[573] = v[223] * v[562] + v[226] * v[565] + v[229] * v[570];
	v[572] = v[222] * v[562] + v[225] * v[565] + v[228] * v[570];
	v[571] = v[221] * v[562] + v[224] * v[565] + v[227] * v[570];
	v[403] = v[401] + v[402];
	v[406] = v[223] * v[393] + v[226] * v[399] + v[229] * v[403];
	v[405] = v[222] * v[393] + v[225] * v[399] + v[228] * v[403];
	v[404] = v[221] * v[393] + v[224] * v[399] + v[227] * v[403];
	v[313] = v[309] + v[312];
	v[2501] = v[15] * v[313];
	v[323] = duiA[0] * v[16] + dduiA[0] * v[17] + v[15] * v[24];
	v[2471] = v[200] * v[323];
	v[324] = duiA[1] * v[16] + dduiA[1] * v[17] + v[15] * v[25];
	v[2478] = v[200] * v[324];
	v[325] = duiA[2] * v[16] + dduiA[2] * v[17] + v[15] * v[26];
	v[2486] = v[200] * v[325];
	v[326] = duiA[3] * v[16] + dduiA[3] * v[17] + v[15] * v[27];
	v[2470] = v[205] * v[326];
	v[327] = duiA[4] * v[16] + dduiA[4] * v[17] + v[15] * v[28];
	v[2477] = v[205] * v[327];
	v[328] = duiA[5] * v[16] + dduiA[5] * v[17] + v[15] * v[29];
	v[2485] = v[205] * v[328];
	v[329] = duiA[6] * v[16] + dduiA[6] * v[17] + v[15] * v[30];
	v[2469] = v[211] * v[329];
	v[330] = duiA[7] * v[16] + dduiA[7] * v[17] + v[15] * v[31];
	v[2476] = v[211] * v[330];
	v[331] = duiA[8] * v[16] + dduiA[8] * v[17] + v[15] * v[32];
	v[2484] = v[211] * v[331];
	v[332] = duiB[0] * v[16] + dduiB[0] * v[17] + v[15] * v[67];
	v[2468] = -(v[302] * v[332]);
	v[333] = duiB[1] * v[16] + dduiB[1] * v[17] + v[15] * v[68];
	v[2475] = -(v[302] * v[333]);
	v[334] = duiB[2] * v[16] + dduiB[2] * v[17] + v[15] * v[69];
	v[2483] = -(v[302] * v[334]);
	v[335] = duiB[3] * v[16] + dduiB[3] * v[17] + v[15] * v[70];
	v[2467] = -(v[307] * v[335]);
	v[336] = duiB[4] * v[16] + dduiB[4] * v[17] + v[15] * v[71];
	v[2474] = -(v[307] * v[336]);
	v[337] = duiB[5] * v[16] + dduiB[5] * v[17] + v[15] * v[72];
	v[2482] = -(v[307] * v[337]);
	v[338] = duiB[6] * v[16] + dduiB[6] * v[17] + v[15] * v[73];
	v[2466] = -(v[313] * v[338]);
	v[1022] = v[2466] + v[2467] + v[2468] + v[2469] + v[2470] + v[2471];
	v[339] = duiB[7] * v[16] + dduiB[7] * v[17] + v[15] * v[74];
	v[2473] = -(v[313] * v[339]);
	v[1023] = v[2473] + v[2474] + v[2475] + v[2476] + v[2477] + v[2478];
	v[340] = duiB[8] * v[16] + dduiB[8] * v[17] + v[15] * v[75];
	v[2481] = -(v[313] * v[340]);
	v[1008] = v[2481] + v[2482] + v[2483] + v[2484] + v[2485] + v[2486];
	v[341] = v[119] * v[200] + v[122] * v[205] + v[125] * v[211] - v[221] * v[302] - v[224] * v[307] - v[227] * v[313];
	v[342] = v[120] * v[200] + v[123] * v[205] + v[126] * v[211] - v[222] * v[302] - v[225] * v[307] - v[228] * v[313];
	v[343] = v[121] * v[200] + v[124] * v[205] + v[127] * v[211] - v[223] * v[302] - v[226] * v[307] - v[229] * v[313];
	v[351] = sqrt((v[341] * v[341]) + (v[342] * v[342]) + (v[343] * v[343]));
	v[1034] = 1e0 / (v[351] * v[351]);
	v[344] = -(v[101] * v[276]) + v[160] * v[52] + v[167] * v[55] + v[174] * v[58] - v[262] * v[95] - v[269] * v[98];
	v[345] = -(v[102] * v[276]) + v[160] * v[53] + v[167] * v[56] + v[174] * v[59] - v[262] * v[96] - v[269] * v[99];
	v[346] = -(v[100] * v[269]) - v[103] * v[276] + v[160] * v[54] + v[167] * v[57] + v[174] * v[60] - v[262] * v[97];
	v[350] = -(*rA) - (*rB) + v[351];
	v[703] = 2e0 * v[13] * v[350];
	v[2472] = v[703] / 2e0;
	if (v[351] > 0.1e-7) { v01 = 1e0 / v[351]; v02 = (-(v01 / v[351])); v03 = (2e0 * v01) / (v[351] * v[351]); }
	else {
		v01 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[351]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[351]) * (0.2399999997e10
			- 0.1199999994e18 * v[351] - 0.3e17 * (v[351] * v[351]))));
		v02 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[351] + 0.6e25 * Power(v[351], 3)
			+ 0.1799999982e26 * (v[351] * v[351]));
		v03 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[351] - 0.3e17 * (v[351] * v[351]));
	};
	v[355] = v03;
	v[356] = v02;
	v[357] = v01;
	v[358] = v[341] * v[357];
	v[2584] = 2e0 * v[358];
	v[751] = v[2472] * v[358];
	v[520] = (v[358] * v[358]);
	v[2552] = -(v[520] * v[703]);
	v[2499] = v[1022] * v[520];
	v[714] = -(v[520] * v[751]);
	v[359] = v[342] * v[357];
	v[2583] = 2e0 * v[359];
	v[2536] = v[323] * v[358] + v[324] * v[359];
	v[2535] = v[326] * v[358] + v[327] * v[359];
	v[2534] = v[329] * v[358] + v[330] * v[359];
	v[2533] = -(v[332] * v[358]) - v[333] * v[359];
	v[2532] = -(v[335] * v[358]) - v[336] * v[359];
	v[2531] = -(v[338] * v[358]) - v[339] * v[359];
	v[2479] = -(v[358] * v[359]);
	v[752] = v[2472] * v[359];
	v[1024] = (v[2466] + v[2467] + v[2468] + v[2469] + v[2470] + v[2471]) * v[358] + (v[2473] + v[2474] + v[2475] + v[2476]
		+ v[2477] + v[2478]) * v[359];
	v[2494] = v[1024] * v[13];
	v[528] = (v[359] * v[359]);
	v[2554] = -(v[528] * v[703]);
	v[2542] = v[520] + v[528];
	v[2498] = v[1023] * v[528];
	v[715] = -(v[528] * v[752]);
	v[527] = v[2479] * v[313];
	v[526] = v[2479] * v[307];
	v[525] = v[2479] * v[302];
	v[524] = -(v[211] * v[2479]);
	v[523] = -(v[205] * v[2479]);
	v[522] = -(v[200] * v[2479]);
	v[360] = v[343] * v[357];
	v[2582] = 2e0 * v[360];
	v[1006] = v[1023] * v[360];
	v[1005] = v[1022] * v[360];
	v[705] = v[2472] * v[360];
	v[2480] = -(v[360] * v[705]);
	v[1135] = v[2480] * v[313];
	v[2739] = v[1135] + v[2554] * v[313];
	v[2738] = v[1135] + v[2552] * v[313];
	v[1131] = v[2480] * v[307];
	v[2737] = v[1131] + v[2554] * v[307];
	v[2736] = v[1131] + v[2552] * v[307];
	v[1127] = v[2480] * v[302];
	v[2735] = v[1127] + v[2554] * v[302];
	v[2734] = v[1127] + v[2552] * v[302];
	v[1123] = -(v[211] * v[2480]);
	v[2733] = v[1123] - v[211] * v[2554];
	v[2732] = v[1123] - v[211] * v[2552];
	v[1119] = -(v[205] * v[2480]);
	v[2731] = v[1119] - v[205] * v[2554];
	v[2730] = v[1119] - v[205] * v[2552];
	v[1115] = -(v[200] * v[2480]);
	v[2729] = v[1115] - v[200] * v[2554];
	v[2728] = v[1115] - v[200] * v[2552];
	v[911] = -(v[2480] * v[358]) - v[714];
	v[909] = -(v[2480] * v[359]) - v[715];
	v[536] = (v[360] * v[360]);
	v[2568] = v[325] * v[536];
	v[2565] = v[328] * v[536];
	v[2562] = v[331] * v[536];
	v[2559] = -(v[334] * v[536]);
	v[2555] = -(v[337] * v[536]);
	v[2550] = -(v[340] * v[536]);
	v[2497] = v[1008] * v[536];
	v[979] = v[2550] + v[2531] * v[360];
	v[960] = v[2555] + v[2532] * v[360];
	v[941] = v[2559] + v[2533] * v[360];
	v[922] = v[2562] + v[2534] * v[360];
	v[900] = v[2565] + v[2535] * v[360];
	v[869] = v[2568] + v[2536] * v[360];
	v[1027] = (v[2481] + v[2482] + v[2483] + v[2484] + v[2485] + v[2486]) * v[360];
	v[706] = v[1027] * v[2472];
	v[361] = sqrt((v[344] * v[344]) + (v[345] * v[345]) + (v[346] * v[346]));
	if (v[361] > 0.1e-7) { v04 = 1e0 / v[361]; v05 = (-(v04 / v[361])); v06 = (2e0 * v04) / (v[361] * v[361]); }
	else {
		v04 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[361]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[361]) * (0.2399999997e10
			- 0.1199999994e18 * v[361] - 0.3e17 * (v[361] * v[361]))));
		v05 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[361] + 0.6e25 * Power(v[361], 3)
			+ 0.1799999982e26 * (v[361] * v[361]));
		v06 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[361] - 0.3e17 * (v[361] * v[361]));
	};
	v[366] = v04;
	v[367] = v[344] * v[366];
	v[368] = v[345] * v[366];
	v[369] = v[346] * v[366];
	v[408] = v[341] * v[374] + v[200] * v[385];
	v[409] = v[342] * v[374] + v[200] * v[386];
	v[410] = v[343] * v[374] + v[200] * v[387];
	v[411] = v[341] * v[380] + v[205] * v[385];
	v[412] = v[342] * v[380] + v[205] * v[386];
	v[413] = v[343] * v[380] + v[205] * v[387];
	v[414] = v[341] * v[384] + v[211] * v[385];
	v[415] = v[342] * v[384] + v[211] * v[386];
	v[416] = v[343] * v[384] + v[211] * v[387];
	v[417] = -(v[302] * v[385]);
	v[418] = -(v[302] * v[386]);
	v[419] = -(v[302] * v[387]);
	v[420] = -(v[307] * v[385]);
	v[421] = -(v[307] * v[386]);
	v[422] = -(v[307] * v[387]);
	v[423] = -(v[313] * v[385]);
	v[424] = -(v[313] * v[386]);
	v[425] = -(v[313] * v[387]);
	v[426] = -(v[200] * v[404]);
	v[427] = -(v[200] * v[405]);
	v[428] = -(v[200] * v[406]);
	v[429] = -(v[205] * v[404]);
	v[430] = -(v[205] * v[405]);
	v[431] = -(v[205] * v[406]);
	v[432] = -(v[211] * v[404]);
	v[433] = -(v[211] * v[405]);
	v[434] = -(v[211] * v[406]);
	v[435] = -(v[341] * v[393]) + v[302] * v[404];
	v[436] = -(v[342] * v[393]) + v[302] * v[405];
	v[437] = -(v[343] * v[393]) + v[302] * v[406];
	v[438] = -(v[341] * v[399]) + v[307] * v[404];
	v[439] = -(v[342] * v[399]) + v[307] * v[405];
	v[440] = -(v[343] * v[399]) + v[307] * v[406];
	v[441] = -(v[341] * v[403]) + v[313] * v[404];
	v[442] = -(v[342] * v[403]) + v[313] * v[405];
	v[443] = -(v[343] * v[403]) + v[313] * v[406];
	v[444] = -(v[18] * v[408]) - v[19] * v[426];
	v[445] = -(v[18] * v[409]) - v[19] * v[427];
	v[446] = -(v[18] * v[410]) - v[19] * v[428];
	v[447] = -(v[18] * v[411]) - v[19] * v[429];
	v[448] = -(v[18] * v[412]) - v[19] * v[430];
	v[449] = -(v[18] * v[413]) - v[19] * v[431];
	v[450] = -(v[18] * v[414]) - v[19] * v[432];
	v[451] = -(v[18] * v[415]) - v[19] * v[433];
	v[452] = -(v[18] * v[416]) - v[19] * v[434];
	v[453] = -(v[18] * v[417]) - v[19] * v[435];
	v[454] = -(v[18] * v[418]) - v[19] * v[436];
	v[455] = -(v[18] * v[419]) - v[19] * v[437];
	v[456] = -(v[18] * v[420]) - v[19] * v[438];
	v[457] = -(v[18] * v[421]) - v[19] * v[439];
	v[458] = -(v[18] * v[422]) - v[19] * v[440];
	v[459] = -(v[18] * v[423]) - v[19] * v[441];
	v[460] = -(v[18] * v[424]) - v[19] * v[442];
	v[461] = -(v[18] * v[425]) - v[19] * v[443];
	v[1302] = -(v[323] * v[444]) - v[324] * v[445] - v[325] * v[446] - v[326] * v[447] - v[327] * v[448] - v[328] * v[449]
		- v[329] * v[450] - v[330] * v[451] - v[331] * v[452] - v[332] * v[453] - v[333] * v[454] - v[334] * v[455] - v[335] * v[456]
		- v[336] * v[457] - v[337] * v[458] - v[338] * v[459] - v[339] * v[460] - v[340] * v[461];
	v[462] = -(v[20] * v[408]) - v[21] * v[426];
	v[1745] = v[1177] * v[444] - v[1193] * v[462];
	v[1726] = v[1176] * v[444] - v[1192] * v[462];
	v[1707] = v[1175] * v[444] - v[1191] * v[462];
	v[1339] = v[444] * v[553] - v[462] * v[571];
	v[1298] = v[444] * v[554] - v[462] * v[572];
	v[1259] = v[444] * v[555] - v[462] * v[573];
	v[463] = -(v[20] * v[409]) - v[21] * v[427];
	v[1746] = v[1177] * v[445] - v[1193] * v[463];
	v[1727] = v[1176] * v[445] - v[1192] * v[463];
	v[1708] = v[1175] * v[445] - v[1191] * v[463];
	v[1337] = v[445] * v[553] - v[463] * v[571];
	v[1296] = v[445] * v[554] - v[463] * v[572];
	v[1257] = v[445] * v[555] - v[463] * v[573];
	v[464] = -(v[20] * v[410]) - v[21] * v[428];
	v[1747] = v[1177] * v[446] - v[1193] * v[464];
	v[1728] = v[1176] * v[446] - v[1192] * v[464];
	v[1709] = v[1175] * v[446] - v[1191] * v[464];
	v[1335] = v[446] * v[553] - v[464] * v[571];
	v[1294] = v[446] * v[554] - v[464] * v[572];
	v[1255] = v[446] * v[555] - v[464] * v[573];
	v[465] = -(v[20] * v[411]) - v[21] * v[429];
	v[1748] = v[1177] * v[447] - v[1193] * v[465];
	v[1729] = v[1176] * v[447] - v[1192] * v[465];
	v[1710] = v[1175] * v[447] - v[1191] * v[465];
	v[1333] = v[447] * v[553] - v[465] * v[571];
	v[1292] = v[447] * v[554] - v[465] * v[572];
	v[1253] = v[447] * v[555] - v[465] * v[573];
	v[466] = -(v[20] * v[412]) - v[21] * v[430];
	v[1749] = v[1177] * v[448] - v[1193] * v[466];
	v[1730] = v[1176] * v[448] - v[1192] * v[466];
	v[1711] = v[1175] * v[448] - v[1191] * v[466];
	v[1331] = v[448] * v[553] - v[466] * v[571];
	v[1290] = v[448] * v[554] - v[466] * v[572];
	v[1251] = v[448] * v[555] - v[466] * v[573];
	v[467] = -(v[20] * v[413]) - v[21] * v[431];
	v[1750] = v[1177] * v[449] - v[1193] * v[467];
	v[1731] = v[1176] * v[449] - v[1192] * v[467];
	v[1712] = v[1175] * v[449] - v[1191] * v[467];
	v[1329] = v[449] * v[553] - v[467] * v[571];
	v[1288] = v[449] * v[554] - v[467] * v[572];
	v[1249] = v[449] * v[555] - v[467] * v[573];
	v[468] = -(v[20] * v[414]) - v[21] * v[432];
	v[1751] = v[1177] * v[450] - v[1193] * v[468];
	v[1732] = v[1176] * v[450] - v[1192] * v[468];
	v[1713] = v[1175] * v[450] - v[1191] * v[468];
	v[1327] = v[450] * v[553] - v[468] * v[571];
	v[1286] = v[450] * v[554] - v[468] * v[572];
	v[1247] = v[450] * v[555] - v[468] * v[573];
	v[469] = -(v[20] * v[415]) - v[21] * v[433];
	v[1752] = v[1177] * v[451] - v[1193] * v[469];
	v[1733] = v[1176] * v[451] - v[1192] * v[469];
	v[1714] = v[1175] * v[451] - v[1191] * v[469];
	v[1325] = v[451] * v[553] - v[469] * v[571];
	v[1284] = v[451] * v[554] - v[469] * v[572];
	v[1245] = v[451] * v[555] - v[469] * v[573];
	v[470] = -(v[20] * v[416]) - v[21] * v[434];
	v[1753] = v[1177] * v[452] - v[1193] * v[470];
	v[1734] = v[1176] * v[452] - v[1192] * v[470];
	v[1715] = v[1175] * v[452] - v[1191] * v[470];
	v[1323] = v[452] * v[553] - v[470] * v[571];
	v[1282] = v[452] * v[554] - v[470] * v[572];
	v[1243] = v[452] * v[555] - v[470] * v[573];
	v[471] = -(v[20] * v[417]) - v[21] * v[435];
	v[1754] = v[1177] * v[453] - v[1193] * v[471];
	v[1735] = v[1176] * v[453] - v[1192] * v[471];
	v[1716] = v[1175] * v[453] - v[1191] * v[471];
	v[1321] = v[453] * v[553] - v[471] * v[571];
	v[1280] = v[453] * v[554] - v[471] * v[572];
	v[1241] = v[453] * v[555] - v[471] * v[573];
	v[472] = -(v[20] * v[418]) - v[21] * v[436];
	v[1755] = v[1177] * v[454] - v[1193] * v[472];
	v[1736] = v[1176] * v[454] - v[1192] * v[472];
	v[1717] = v[1175] * v[454] - v[1191] * v[472];
	v[1319] = v[454] * v[553] - v[472] * v[571];
	v[1278] = v[454] * v[554] - v[472] * v[572];
	v[1239] = v[454] * v[555] - v[472] * v[573];
	v[473] = -(v[20] * v[419]) - v[21] * v[437];
	v[1756] = v[1177] * v[455] - v[1193] * v[473];
	v[1737] = v[1176] * v[455] - v[1192] * v[473];
	v[1718] = v[1175] * v[455] - v[1191] * v[473];
	v[1317] = v[455] * v[553] - v[473] * v[571];
	v[1276] = v[455] * v[554] - v[473] * v[572];
	v[1237] = v[455] * v[555] - v[473] * v[573];
	v[474] = -(v[20] * v[420]) - v[21] * v[438];
	v[1757] = v[1177] * v[456] - v[1193] * v[474];
	v[1738] = v[1176] * v[456] - v[1192] * v[474];
	v[1719] = v[1175] * v[456] - v[1191] * v[474];
	v[1315] = v[456] * v[553] - v[474] * v[571];
	v[1274] = v[456] * v[554] - v[474] * v[572];
	v[1235] = v[456] * v[555] - v[474] * v[573];
	v[475] = -(v[20] * v[421]) - v[21] * v[439];
	v[1758] = v[1177] * v[457] - v[1193] * v[475];
	v[1739] = v[1176] * v[457] - v[1192] * v[475];
	v[1720] = v[1175] * v[457] - v[1191] * v[475];
	v[1313] = v[457] * v[553] - v[475] * v[571];
	v[1272] = v[457] * v[554] - v[475] * v[572];
	v[1233] = v[457] * v[555] - v[475] * v[573];
	v[476] = -(v[20] * v[422]) - v[21] * v[440];
	v[1759] = v[1177] * v[458] - v[1193] * v[476];
	v[1740] = v[1176] * v[458] - v[1192] * v[476];
	v[1721] = v[1175] * v[458] - v[1191] * v[476];
	v[1311] = v[458] * v[553] - v[476] * v[571];
	v[1270] = v[458] * v[554] - v[476] * v[572];
	v[1231] = v[458] * v[555] - v[476] * v[573];
	v[477] = -(v[20] * v[423]) - v[21] * v[441];
	v[1760] = v[1177] * v[459] - v[1193] * v[477];
	v[1741] = v[1176] * v[459] - v[1192] * v[477];
	v[1722] = v[1175] * v[459] - v[1191] * v[477];
	v[1309] = v[459] * v[553] - v[477] * v[571];
	v[1268] = v[459] * v[554] - v[477] * v[572];
	v[1229] = v[459] * v[555] - v[477] * v[573];
	v[478] = -(v[20] * v[424]) - v[21] * v[442];
	v[1761] = v[1177] * v[460] - v[1193] * v[478];
	v[1742] = v[1176] * v[460] - v[1192] * v[478];
	v[1723] = v[1175] * v[460] - v[1191] * v[478];
	v[1307] = v[460] * v[553] - v[478] * v[571];
	v[1266] = v[460] * v[554] - v[478] * v[572];
	v[1227] = v[460] * v[555] - v[478] * v[573];
	v[479] = -(v[20] * v[425]) - v[21] * v[443];
	v[1762] = v[1177] * v[461] - v[1193] * v[479];
	v[1743] = v[1176] * v[461] - v[1192] * v[479];
	v[1724] = v[1175] * v[461] - v[1191] * v[479];
	v[1305] = v[461] * v[553] - v[479] * v[571];
	v[1300] = v[323] * v[462] + v[324] * v[463] + v[325] * v[464] + v[326] * v[465] + v[327] * v[466] + v[328] * v[467]
		+ v[329] * v[468] + v[330] * v[469] + v[331] * v[470] + v[332] * v[471] + v[333] * v[472] + v[334] * v[473] + v[335] * v[474]
		+ v[336] * v[475] + v[337] * v[476] + v[338] * v[477] + v[339] * v[478] + v[340] * v[479];
	v[1264] = v[461] * v[554] - v[479] * v[572];
	v[1225] = v[461] * v[555] - v[479] * v[573];
	b480 = sqrt(Power(v[359] * v[367] - v[358] * v[368], 2) + Power(-(v[360] * v[367]) + v[358] * v[369], 2) + Power
	(v[360] * v[368] - v[359] * v[369], 2)) > 0.1e-7;
	if (b480) {
		v[482] = v[360] * v[368] - v[359] * v[369];
		v[483] = -(v[360] * v[367]) + v[358] * v[369];
		v[484] = v[359] * v[367] - v[358] * v[368];
		v[485] = sqrt((v[482] * v[482]) + (v[483] * v[483]) + (v[484] * v[484]));
		v[1626] = 1e0 / (v[485] * v[485]);
		v[1436] = v[485];
		v[1637] = 1e0 - (v[1436] * v[1436]);
		v[2773] = 1e0 / Power(v[1637], 0.15e1);
		v[1632] = 1e0 / sqrt(v[1637]);
		v[1435] = asin(v[1436]) / 2e0;
		v[1631] = 1e0 / Power(cos(v[1435]), 2);
		v[2605] = v[1631] * v[1632];
		v[487] = 2e0 * tan(v[1435]);
		if (v[485] > 0.1e-7) { v07 = 1e0 / v[485]; v08 = (-(v07 / v[485])); v09 = (2e0 * v07) / (v[485] * v[485]); }
		else {
			v07 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[485]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[485]) *
				(0.2399999997e10 - 0.1199999994e18 * v[485] - 0.3e17 * (v[485] * v[485]))));
			v08 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[485] + 0.6e25 * Power(v[485], 3)
				+ 0.1799999982e26 * (v[485] * v[485]));
			v09 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[485] - 0.3e17 * (v[485] * v[485]));
		};
		v[491] = v09;
		v[492] = v08;
		v[493] = v07;
		v[2772] = v[487] * v[492] + v[2605] * v[493];
		v[2487] = v[487] * v[493];
		v[494] = v[2487] * v[482];
		v[2643] = 2e0 * v[494];
		v[2579] = v[494] / 2e0;
		v[505] = (v[494] * v[494]);
		v[495] = v[2487] * v[483];
		v[2488] = v[495] / 2e0;
		v[503] = v[2488] * v[494];
		v[498] = (v[495] * v[495]);
		v[1390] = -v[498] - v[505];
		v[496] = v[2487] * v[484];
		v[2640] = 2e0 * v[496];
		v[1420] = -v[496] + v[503];
		v[1411] = v[496] + v[503];
		v[510] = v[2488] * v[496];
		v[1402] = -v[494] + v[510];
		v[1394] = v[494] + v[510];
		v[508] = v[2579] * v[496];
		v[1415] = v[495] + v[508];
		v[1398] = -v[495] + v[508];
		v[499] = (v[496] * v[496]);
		v[1430] = 4e0 + v[498] + v[499] + v[505];
		v[2774] = 1e0 / Power(v[1430], 3);
		v[2639] = -4e0 / (v[1430] * v[1430]);
		v[1425] = -v[498] - v[499];
		v[1407] = -v[499] - v[505];
		v[497] = 4e0 / v[1430];
		v[2489] = v[497] / 2e0;
		v[500] = 1e0 + v[1425] * v[2489];
		v[501] = v[1420] * v[497];
		v[502] = v[1415] * v[497];
		v[504] = v[1411] * v[497];
		v[506] = 1e0 + v[1407] * v[2489];
		v[507] = v[1402] * v[497];
		v[509] = v[1398] * v[497];
		v[511] = v[1394] * v[497];
		v[512] = 1e0 + v[1390] * v[2489];
	}
	else {
		v[500] = 1e0;
		v[501] = 0e0;
		v[502] = 0e0;
		v[504] = 0e0;
		v[506] = 1e0;
		v[507] = 0e0;
		v[509] = 0e0;
		v[511] = 0e0;
		v[512] = 1e0;
	};
	if (b23) {
		v[1365] = 1e0 - v[536];
		v[1363] = 1e0 - v[528];
		v[1361] = 1e0 - v[520];
		v[517] = v[121] * v[160] + v[124] * v[167] + v[127] * v[174] - v[223] * v[262] - v[226] * v[269] - v[229] * v[276]
			+ v[5] * v[509] + v[511] * v[6] + v[512] * v[7];
		v[2490] = v[360] * v[517];
		v[516] = v[120] * v[160] + v[123] * v[167] + v[126] * v[174] - v[222] * v[262] - v[225] * v[269] - v[228] * v[276]
			+ v[5] * v[504] + v[506] * v[6] + v[507] * v[7];
		v[2492] = v[359] * v[516];
		v[2578] = v[2490] + v[2492];
		v[515] = v[119] * v[160] + v[122] * v[167] + v[125] * v[174] - v[221] * v[262] - v[224] * v[269] - v[227] * v[276]
			+ v[5] * v[500] + v[501] * v[6] + v[502] * v[7];
		v[2491] = -(v[358] * v[515]);
		v[2577] = -v[2490] + v[2491];
		v[2576] = v[2491] - v[2492];
		v[514] = -(v[2578] * v[358]) + v[1361] * v[515];
		v[518] = v[2577] * v[359] + v[1363] * v[516];
		v[519] = v[2576] * v[360] + v[1365] * v[517];
	}
	else {
		v[514] = 0e0;
		v[518] = 0e0;
		v[519] = 0e0;
	};
	v[577] = fabs(v[350]);
	v[2493] = (*epsn) / v[577];
	v[1038] = Power(v[577], v[9]);
	v[1037] = v[2493] * v[350];
	v[579] = v[1037] * v[1038];
	v[2649] = v[1038] * v[2493] + _copysign(1.e0, v[350]) * (-(v[579] / v[577]) + v[1037] * v[9] * Power(v[577], -1e0
		+ v[9]));
	v[1920] = v[2494] + v[579];
	v[1031] = v[1027] * v[13] + v[579];
	v[582] = v[13] * (v[2499] + v[1027] * v[358] + v[324] * v[522] + v[327] * v[523] + v[330] * v[524] + v[333] * v[525]
		+ v[336] * v[526] + v[339] * v[527]);
	v[583] = v[13] * (v[2498] + v[1027] * v[359] + v[323] * v[522] + v[326] * v[523] + v[329] * v[524] + v[332] * v[525]
		+ v[335] * v[526] + v[338] * v[527]);
	v[584] = v[13] * (v[2497] + v[1024] * v[360]);
	v[585] = v[358] * v[579] + v[582];
	v[586] = v[359] * v[579] + v[583];
	v[587] = v[360] * v[579] + v[584];
	v[2547] = v[2494] * v[360] + v[587];
	v[1791] = (v[585] * v[585]) + (v[586] * v[586]) + (v[587] * v[587]);
	v[588] = v[10] * v[514];
	v[589] = v[10] * v[518];
	v[590] = v[10] * v[519];
	v[594] = -(v[14] * (v[1339] * v[323] + v[1337] * v[324] + v[1335] * v[325] + v[1333] * v[326] + v[1331] * v[327]
		+ v[1329] * v[328] + v[1327] * v[329] + v[1325] * v[330] + v[1323] * v[331] + v[1321] * v[332] + v[1319] * v[333]
		+ v[1317] * v[334] + v[1315] * v[335] + v[1313] * v[336] + v[1311] * v[337] + v[1309] * v[338] + v[1307] * v[339]
		+ v[1305] * v[340])) + v[588];
	v[595] = -(v[14] * (v[1298] * v[323] + v[1296] * v[324] + v[1294] * v[325] + v[1292] * v[326] + v[1290] * v[327]
		+ v[1288] * v[328] + v[1286] * v[329] + v[1284] * v[330] + v[1282] * v[331] + v[1280] * v[332] + v[1278] * v[333]
		+ v[1276] * v[334] + v[1274] * v[335] + v[1272] * v[336] + v[1270] * v[337] + v[1268] * v[338] + v[1266] * v[339]
		+ v[1264] * v[340])) + v[589];
	v[596] = -(v[14] * (v[1259] * v[323] + v[1257] * v[324] + v[1255] * v[325] + v[1253] * v[326] + v[1251] * v[327]
		+ v[1249] * v[328] + v[1247] * v[329] + v[1245] * v[330] + v[1243] * v[331] + v[1241] * v[332] + v[1239] * v[333]
		+ v[1237] * v[334] + v[1235] * v[335] + v[1233] * v[336] + v[1231] * v[337] + v[1229] * v[338] + v[1227] * v[339]
		+ v[1225] * v[340])) + v[590];
	v[1787] = (v[594] * v[594]) + (v[595] * v[595]) + (v[596] * v[596]);
	if (b22) {
		b598 = sqrt((v[594] * v[594]) + (v[595] * v[595]) + (v[596] * v[596])) <= v[11] * sqrt((v[585] * v[585]) +
			(v[586] * v[586]) + (v[587] * v[587]));
		if (b598) {
			v[600] = v[594];
			v[601] = v[595];
			v[602] = v[596];
			v[603] = 1e0;
		}
		else {
			v[2495] = v[12] * sqrt(v[1791]);
			v[604] = sqrt(v[1787]);
			if (v[604] > 0.1e-5) { v010 = 1e0 / v[604]; v011 = (-(v010 / v[604])); v012 = (2e0 * v010) / (v[604] * v[604]); }
			else {
				v010 = (24000000e0 - (-1e0 + 1000000e0 * v[604]) * (71999994e0 - 0.71999982e14 * v[604] + 0.6e19 * Power(v[604], 3)
					+ 0.23999982e20 * (v[604] * v[604]))) / 24e0;
				v011 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[604] + 0.6e19 * Power(v[604], 3) + 0.17999982e20 *
					(v[604] * v[604]));
				v012 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[604] - 0.3e13 * (v[604] * v[604]));
			};
			v[608] = v011;
			v[609] = v010;
			v[610] = v[594] * v[609];
			v[611] = v[595] * v[609];
			v[612] = v[596] * v[609];
			v[600] = v[2495] * v[610];
			v[601] = v[2495] * v[611];
			v[602] = v[2495] * v[612];
			v[603] = 0e0;
		};
		if (sqrt((v[588] * v[588]) + (v[589] * v[589]) + (v[590] * v[590])) > v[11] * sqrt((v[585] * v[585]) + (v[586] * v[586]
			) + (v[587] * v[587]))) {
			if (v[10] > 0.1e-5) { v013 = 1e0 / v[10]; v014 = (-(v013 / v[10])); v015 = (2e0 * v013) / (v[10] * v[10]); }
			else {
				v013 = (24000000e0 - (-1e0 + 1000000e0 * v[10]) * (71999994e0 - 0.71999982e14 * v[10] + 0.6e19 * Power(v[10], 3)
					+ 0.23999982e20 * (v[10] * v[10]))) / 24e0;
				v014 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[10] + 0.6e19 * Power(v[10], 3) + 0.17999982e20 *
					(v[10] * v[10]));
				v015 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[10] - 0.3e13 * (v[10] * v[10]));
			};
			v[622] = sqrt((v[588] * v[588]) + (v[589] * v[589]) + (v[590] * v[590]));
			if (v[622] > 0.1e-5) { v016 = 1e0 / v[622]; v017 = (-(v016 / v[622])); v018 = (2e0 * v016) / (v[622] * v[622]); }
			else {
				v016 = (24000000e0 - (-1e0 + 1000000e0 * v[622]) * (71999994e0 - 0.71999982e14 * v[622] + 0.6e19 * Power(v[622], 3)
					+ 0.23999982e20 * (v[622] * v[622]))) / 24e0;
				v017 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[622] + 0.6e19 * Power(v[622], 3) + 0.17999982e20 *
					(v[622] * v[622]));
				v018 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[622] - 0.3e13 * (v[622] * v[622]));
			};
			v[629] = -(v013 * v016 * v[12] * sqrt(v[1791]));
			v[628] = v[514] + v[588] * v[629];
			v[630] = v[518] + v[589] * v[629];
			v[631] = v[519] + v[590] * v[629];
		}
		else {
			v[628] = 0e0;
			v[630] = 0e0;
			v[631] = 0e0;
		};
	}
	else {
		b632 = sqrt((v[594] * v[594]) + (v[595] * v[595]) + (v[596] * v[596])) <= v[12] * sqrt((v[585] * v[585]) +
			(v[586] * v[586]) + (v[587] * v[587]));
		if (b632) {
			v[600] = v[594];
			v[601] = v[595];
			v[602] = v[596];
			v[603] = 1e0;
		}
		else {
			v[643] = sqrt(v[1791]);
			v[2496] = v[12] * v[643];
			v[634] = sqrt(v[1787]);
			if (v[634] > 0.1e-5) { v019 = 1e0 / v[634]; v020 = (-(v019 / v[634])); v021 = (2e0 * v019) / (v[634] * v[634]); }
			else {
				v019 = (24000000e0 - (-1e0 + 1000000e0 * v[634]) * (71999994e0 - 0.71999982e14 * v[634] + 0.6e19 * Power(v[634], 3)
					+ 0.23999982e20 * (v[634] * v[634]))) / 24e0;
				v020 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[634] + 0.6e19 * Power(v[634], 3) + 0.17999982e20 *
					(v[634] * v[634]));
				v021 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[634] - 0.3e13 * (v[634] * v[634]));
			};
			v[638] = v020;
			v[639] = v019;
			v[640] = v[594] * v[639];
			v[641] = v[595] * v[639];
			v[642] = v[596] * v[639];
			v[600] = v[2496] * v[640];
			v[601] = v[2496] * v[641];
			v[602] = v[2496] * v[642];
			v[603] = 0e0;
		};
		if (sqrt((v[588] * v[588]) + (v[589] * v[589]) + (v[590] * v[590])) > v[12] * sqrt((v[585] * v[585]) + (v[586] * v[586]
			) + (v[587] * v[587]))) {
			if (v[10] > 0.1e-5) { v022 = 1e0 / v[10]; v023 = (-(v022 / v[10])); v024 = (2e0 * v022) / (v[10] * v[10]); }
			else {
				v022 = (24000000e0 - (-1e0 + 1000000e0 * v[10]) * (71999994e0 - 0.71999982e14 * v[10] + 0.6e19 * Power(v[10], 3)
					+ 0.23999982e20 * (v[10] * v[10]))) / 24e0;
				v023 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[10] + 0.6e19 * Power(v[10], 3) + 0.17999982e20 *
					(v[10] * v[10]));
				v024 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[10] - 0.3e13 * (v[10] * v[10]));
			};
			v[652] = sqrt((v[588] * v[588]) + (v[589] * v[589]) + (v[590] * v[590]));
			if (v[652] > 0.1e-5) { v025 = 1e0 / v[652]; v026 = (-(v025 / v[652])); v027 = (2e0 * v025) / (v[652] * v[652]); }
			else {
				v025 = (24000000e0 - (-1e0 + 1000000e0 * v[652]) * (71999994e0 - 0.71999982e14 * v[652] + 0.6e19 * Power(v[652], 3)
					+ 0.23999982e20 * (v[652] * v[652]))) / 24e0;
				v026 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[652] + 0.6e19 * Power(v[652], 3) + 0.17999982e20 *
					(v[652] * v[652]));
				v027 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[652] - 0.3e13 * (v[652] * v[652]));
			};
			v[658] = -(v022 * v025 * v[12] * sqrt(v[1791]));
			v[628] = v[514] + v[588] * v[658];
			v[630] = v[518] + v[589] * v[658];
			v[631] = v[519] + v[590] * v[658];
		}
		else {
			v[628] = 0e0;
			v[630] = 0e0;
			v[631] = 0e0;
		};
	};
	fn[0] = v[585];
	fn[1] = v[586];
	fn[2] = v[587];
	ft[0] = v[600];
	ft[1] = v[601];
	ft[2] = v[602];
	(*stickupdated) = v[603];
	gtpupdated[0] = v[514] - v[628];
	gtpupdated[1] = v[518] - v[630];
	gtpupdated[2] = v[519] - v[631];
	v[669] = v[2472] * v[2542];
	v[2697] = v[669] + v[536] * v[703];
	v[1086] = v[325] * v[669];
	v[2721] = v[1086] + v[2568] * v[703];
	v[2703] = v[1086] + v[2536] * v[705];
	v[1080] = v[328] * v[669];
	v[2718] = v[1080] + v[2565] * v[703];
	v[2702] = v[1080] + v[2535] * v[705];
	v[1074] = v[331] * v[669];
	v[2715] = v[1074] + v[2562] * v[703];
	v[2701] = v[1074] + v[2534] * v[705];
	v[1068] = -(v[334] * v[669]);
	v[2712] = v[1068] + v[2559] * v[703];
	v[2700] = v[1068] + v[2533] * v[705];
	v[1059] = -(v[337] * v[669]);
	v[2709] = v[1059] + v[2555] * v[703];
	v[2699] = v[1059] + v[2532] * v[705];
	v[1050] = -(v[340] * v[669]);
	v[2704] = v[1050] + v[2550] * v[703];
	v[2698] = v[1050] + v[2531] * v[705];
	v[907] = v[360] * v[669] + v[536] * v[705];
	v[670] = v[324] * v[751] + v[323] * v[752];
	v[2567] = -(v[2480] * v[324]) + v[358] * v[670];
	v[2720] = v[2567] - v[2554] * v[324];
	v[1082] = v[359] * v[670];
	v[2566] = v[1082] - v[2480] * v[323];
	v[2719] = v[2566] - v[2552] * v[323];
	v[673] = v[327] * v[751] + v[326] * v[752];
	v[2564] = -(v[2480] * v[327]) + v[358] * v[673];
	v[2717] = v[2564] - v[2554] * v[327];
	v[1076] = v[359] * v[673];
	v[2563] = v[1076] - v[2480] * v[326];
	v[2716] = v[2563] - v[2552] * v[326];
	v[674] = v[330] * v[751] + v[329] * v[752];
	v[2561] = -(v[2480] * v[330]) + v[358] * v[674];
	v[2714] = v[2561] - v[2554] * v[330];
	v[1070] = v[359] * v[674];
	v[2560] = v[1070] - v[2480] * v[329];
	v[2713] = v[2560] - v[2552] * v[329];
	v[675] = v[333] * v[751] + v[332] * v[752];
	v[2558] = v[2480] * v[333] - v[358] * v[675];
	v[2711] = v[2558] + v[2554] * v[333];
	v[1064] = -(v[359] * v[675]);
	v[2557] = v[1064] + v[2480] * v[332];
	v[2710] = v[2557] + v[2552] * v[332];
	v[676] = v[336] * v[751] + v[335] * v[752];
	v[2553] = v[2480] * v[336] - v[358] * v[676];
	v[2708] = v[2553] + v[2554] * v[336];
	v[1055] = -(v[359] * v[676]);
	v[2551] = v[1055] + v[2480] * v[335];
	v[2707] = v[2551] + v[2552] * v[335];
	v[677] = v[339] * v[751] + v[338] * v[752];
	v[2549] = v[2480] * v[339] - v[358] * v[677];
	v[2706] = v[2549] + v[2554] * v[339];
	v[1046] = -(v[359] * v[677]);
	v[2548] = v[1046] + v[2480] * v[338];
	v[2705] = v[2548] + v[2552] * v[338];
	v[1030] = v[200] * v[670] + v[205] * v[673] + v[211] * v[674] - v[302] * v[675] - v[307] * v[676] - v[313] * v[677];
	v[702] = v[2547] * v[350] + v[1008] * v[669] + v[2497] * v[703];
	v[704] = v[1030] * v[358] + v[350] * v[586] + v[2498] * v[703] + v[1006] * v[705] + v[359] * v[706];
	v[707] = v[1030] * v[359] + v[350] * v[585] + v[2499] * v[703] + v[1005] * v[705] + v[358] * v[706];
	v[1035] = v[343] * v[702] + v[342] * v[704] + v[341] * v[707];
	v[708] = v[1035] * v[356] + v[358] * v[585] + v[359] * v[586] + v[360] * v[587];
	v[2500] = v[708] / v[351];
	v[709] = v[2500] * v[343] + v[357] * v[702];
	v[1052] = -v[709] - v[15] * v[907];
	v[711] = v[2500] * v[342] + v[357] * v[704];
	v[1053] = -v[711] - v[15] * v[909];
	v[712] = v[2500] * v[341] + v[357] * v[707];
	v[1054] = -v[712] - v[15] * v[911];
	v[713] = v[1046] * v[358] + v[1050] * v[360] - v[229] * v[709] - v[228] * v[711] - v[227] * v[712] + v[338] * v[714]
		+ v[339] * v[715] + v[705] * v[979];
	v[716] = v[1055] * v[358] + v[1059] * v[360] - v[226] * v[709] - v[225] * v[711] - v[224] * v[712] + v[335] * v[714]
		+ v[336] * v[715] + v[705] * v[960];
	v[717] = v[1064] * v[358] + v[1068] * v[360] - v[223] * v[709] - v[222] * v[711] - v[221] * v[712] + v[332] * v[714]
		+ v[333] * v[715] + v[705] * v[941];
	v[718] = v[1070] * v[358] + v[1074] * v[360] + v[127] * v[709] + v[126] * v[711] + v[125] * v[712] - v[329] * v[714]
		- v[330] * v[715] + v[705] * v[922];
	v[719] = v[1076] * v[358] + v[1080] * v[360] + v[124] * v[709] + v[123] * v[711] + v[122] * v[712] - v[326] * v[714]
		- v[327] * v[715] + v[705] * v[900];
	v[720] = v[1082] * v[358] + v[1086] * v[360] + v[121] * v[709] + v[120] * v[711] + v[119] * v[712] - v[323] * v[714]
		- v[324] * v[715] + v[705] * v[869];
	v[721] = 0e0;
	b722 = b273;
	if (b722) {
		v[721] = v[568] * v[713];
	}
	else {
	};
	v[723] = 0e0;
	b724 = b266;
	if (b724) {
		v[723] = v[2465] * v[713];
		v[721] = -(v[400] * v[713]) + v[721];
	}
	else {
	};
	b726 = b266;
	if (b726) {
		v[723] = v[2464] * v[716] + v[723];
		v[721] = v[400] * v[716] + v[721];
	}
	else {
	};
	v[727] = 0e0;
	b728 = b259;
	if (b728) {
		v[727] = v[2463] * v[716];
		v[721] = -(v[394] * v[716]) + v[721];
	}
	else {
	};
	b729 = b259;
	if (b729) {
		v[727] = v[2462] * v[717] + v[727];
		v[721] = v[394] * v[717] + v[721];
	}
	else {
	};
	b731 = b256;
	if (b731) {
		v[721] = v[558] * v[717] + v[721];
	}
	else {
	};
	v[732] = v[723];
	b733 = b239;
	if (b733) {
		v[721] = v[721] + v[244] * v[732];
	}
	else {
	};
	v[734] = v[727];
	b735 = b239;
	if (b735) {
		v[721] = v[721] - v[244] * v[734];
	}
	else {
	};
	v[736] = 0e0;
	b737 = b171;
	if (b737) {
		v[736] = v[550] * v[718];
	}
	else {
	};
	v[738] = 0e0;
	b739 = b164;
	if (b739) {
		v[738] = v[2459] * v[718];
		v[736] = -(v[381] * v[718]) + v[736];
	}
	else {
	};
	b741 = b164;
	if (b741) {
		v[738] = v[2458] * v[719] + v[738];
		v[736] = v[381] * v[719] + v[736];
	}
	else {
	};
	v[742] = 0e0;
	b743 = b157;
	if (b743) {
		v[742] = v[2457] * v[719];
		v[736] = -(v[375] * v[719]) + v[736];
	}
	else {
	};
	b744 = b157;
	if (b744) {
		v[742] = v[2456] * v[720] + v[742];
		v[736] = v[375] * v[720] + v[736];
	}
	else {
	};
	b746 = b154;
	if (b746) {
		v[736] = v[540] * v[720] + v[736];
	}
	else {
	};
	v[747] = v[738];
	b748 = b137;
	if (b748) {
		v[736] = v[736] + v[142] * v[747];
	}
	else {
	};
	v[749] = v[742];
	b750 = b137;
	if (b750) {
		v[736] = v[736] - v[142] * v[749];
	}
	else {
	};
	v[4056] = v[200] * v[712] + v[462] * v[721] + v[444] * v[736] + v[15] * (v[522] * v[752] + v[200] * v[911]);
	v[4057] = v[200] * v[711] + v[463] * v[721] + v[445] * v[736] + v[15] * (v[522] * v[751] + v[200] * v[909]);
	v[4058] = -(v[1052] * v[200]) + v[464] * v[721] + v[446] * v[736];
	v[4059] = v[205] * v[712] + v[465] * v[721] + v[447] * v[736] + v[15] * (v[523] * v[752] + v[205] * v[911]);
	v[4060] = v[205] * v[711] + v[466] * v[721] + v[448] * v[736] + v[15] * (v[523] * v[751] + v[205] * v[909]);
	v[4061] = -(v[1052] * v[205]) + v[467] * v[721] + v[449] * v[736];
	v[4062] = v[211] * v[712] + v[468] * v[721] + v[450] * v[736] + v[15] * (v[524] * v[752] + v[211] * v[911]);
	v[4063] = v[211] * v[711] + v[469] * v[721] + v[451] * v[736] + v[15] * (v[524] * v[751] + v[211] * v[909]);
	v[4064] = -(v[1052] * v[211]) + v[470] * v[721] + v[452] * v[736];
	v[4065] = -(v[302] * v[712]) + v[471] * v[721] + v[453] * v[736] + v[15] * (v[525] * v[752] - v[302] * v[911]);
	v[4066] = -(v[302] * v[711]) + v[472] * v[721] + v[454] * v[736] + v[15] * (v[525] * v[751] - v[302] * v[909]);
	v[4067] = v[1052] * v[302] + v[473] * v[721] + v[455] * v[736];
	v[4068] = -(v[307] * v[712]) + v[474] * v[721] + v[456] * v[736] + v[15] * (v[526] * v[752] - v[307] * v[911]);
	v[4069] = -(v[307] * v[711]) + v[475] * v[721] + v[457] * v[736] + v[15] * (v[526] * v[751] - v[307] * v[909]);
	v[4070] = v[1052] * v[307] + v[476] * v[721] + v[458] * v[736];
	v[4071] = -(v[313] * v[712]) + v[477] * v[721] + v[459] * v[736] + v[15] * (v[527] * v[752] - v[313] * v[911]);
	v[4072] = -(v[313] * v[711]) + v[478] * v[721] + v[460] * v[736] + v[15] * (v[527] * v[751] - v[313] * v[909]);
	v[4073] = v[1052] * v[313] + v[479] * v[721] + v[461] * v[736];
	for (i667 = 1; i667 <= 18; i667++) {
		i2524 = (i667 == 18 ? 1 : 0);
		i2523 = (i667 == 17 ? 1 : 0);
		i2522 = (i667 == 16 ? 1 : 0);
		i2521 = (i667 == 15 ? 1 : 0);
		i2520 = (i667 == 14 ? 1 : 0);
		i2519 = (i667 == 13 ? 1 : 0);
		i2518 = (i667 == 12 ? 1 : 0);
		i2517 = (i667 == 11 ? 1 : 0);
		i2516 = (i667 == 10 ? 1 : 0);
		i2515 = (i667 == 9 ? 1 : 0);
		i2514 = (i667 == 8 ? 1 : 0);
		i2513 = (i667 == 7 ? 1 : 0);
		i2512 = (i667 == 6 ? 1 : 0);
		i2511 = (i667 == 5 ? 1 : 0);
		i2510 = (i667 == 4 ? 1 : 0);
		i2509 = (i667 == 3 ? 1 : 0);
		i2508 = (i667 == 2 ? 1 : 0);
		i2507 = (i667 == 1 ? 1 : 0);
		v[872] = i2523 * v[2501];
		v[880] = i2522 * v[2501];
		v[873] = i2520 * v[2502];
		v[881] = i2519 * v[2502];
		v[874] = i2517 * v[2503];
		v[882] = i2516 * v[2503];
		v[875] = i2514 * v[2504];
		v[883] = i2513 * v[2504];
		v[876] = i2511 * v[2505];
		v[884] = i2510 * v[2505];
		v[2527] = -(i2524 * v[2501]) - i2521 * v[2502] - i2518 * v[2503] + i2515 * v[2504] + i2512 * v[2505] + i2509 * v[2506];
		v[877] = i2508 * v[2506];
		v[2537] = -v[872] - v[873] - v[874] + v[875] + v[876] + v[877];
		v[885] = i2507 * v[2506];
		v[2539] = -v[880] - v[881] - v[882] + v[883] + v[884] + v[885];
		v[813] = i2507 * v[444] + i2508 * v[445] + i2509 * v[446] + i2510 * v[447] + i2511 * v[448] + i2512 * v[449] + i2513 * v[450]
			+ i2514 * v[451] + i2515 * v[452] + i2516 * v[453] + i2517 * v[454] + i2518 * v[455] + i2519 * v[456] + i2520 * v[457]
			+ i2521 * v[458] + i2522 * v[459] + i2523 * v[460] + i2524 * v[461];
		v[814] = i2507 * v[462] + i2508 * v[463] + i2509 * v[464] + i2510 * v[465] + i2511 * v[466] + i2512 * v[467] + i2513 * v[468]
			+ i2514 * v[469] + i2515 * v[470] + i2516 * v[471] + i2517 * v[472] + i2518 * v[473] + i2519 * v[474] + i2520 * v[475]
			+ i2521 * v[476] + i2522 * v[477] + i2523 * v[478] + i2524 * v[479];
		b815 = b137;
		if (b815) {
			v[816] = -(v[142] * v[813]);
		}
		else {
			v[816] = 0e0;
		};
		v[817] = v[816];
		b818 = b137;
		if (b818) {
			v[819] = v[142] * v[813];
		}
		else {
			v[819] = 0e0;
		};
		v[820] = v[819];
		v[821] = 0e0;
		b822 = b154;
		if (b822) {
			v[821] = v[540] * v[813];
		}
		else {
		};
		v[823] = 0e0;
		v[824] = 0e0;
		b825 = b157;
		if (b825) {
			v[821] = v[375] * v[813] + v[821];
			v[823] = v[720] * v[813];
			v[821] = v[2456] * v[817] + v[821];
			v[824] = v[162] * v[720] * v[817];
		}
		else {
		};
		v[826] = 0e0;
		v[827] = 0e0;
		b828 = b157;
		if (b828) {
			v[827] = -(v[375] * v[813]);
			v[823] = -(v[719] * v[813]) + v[823];
			v[827] = v[2457] * v[817] + v[827];
			v[826] = v[162] * v[719] * v[817];
			v[817] = 0e0;
		}
		else {
		};
		v[829] = 0e0;
		v[830] = 0e0;
		b831 = b164;
		if (b831) {
			v[827] = v[381] * v[813] + v[827];
			v[830] = v[719] * v[813];
			v[827] = v[2458] * v[820] + v[827];
			v[829] = v[169] * v[719] * v[820];
		}
		else {
		};
		v[832] = 0e0;
		v[833] = 0e0;
		b834 = b164;
		if (b834) {
			v[832] = -(v[381] * v[813]);
			v[830] = -(v[718] * v[813]) + v[830];
			v[832] = v[2459] * v[820] + v[832];
			v[833] = v[169] * v[718] * v[820];
			v[820] = 0e0;
		}
		else {
		};
		b835 = b171;
		if (b835) {
			v[832] = v[550] * v[813] + v[832];
			v[813] = 0e0;
		}
		else {
		};
		b836 = b239;
		if (b836) {
			v[837] = -(v[244] * v[814]);
		}
		else {
			v[837] = 0e0;
		};
		v[838] = v[837];
		b839 = b239;
		if (b839) {
			v[840] = v[244] * v[814];
		}
		else {
			v[840] = 0e0;
		};
		v[841] = v[840];
		v[842] = 0e0;
		b843 = b256;
		if (b843) {
			v[842] = v[558] * v[814];
		}
		else {
		};
		v[844] = 0e0;
		v[845] = 0e0;
		b846 = b259;
		if (b846) {
			v[842] = v[394] * v[814] + v[842];
			v[844] = v[717] * v[814];
			v[842] = v[2462] * v[838] + v[842];
			v[845] = v[264] * v[717] * v[838];
		}
		else {
		};
		v[847] = 0e0;
		v[848] = 0e0;
		b849 = b259;
		if (b849) {
			v[848] = -(v[394] * v[814]);
			v[844] = -(v[716] * v[814]) + v[844];
			v[848] = v[2463] * v[838] + v[848];
			v[847] = v[264] * v[716] * v[838];
			v[838] = 0e0;
		}
		else {
		};
		v[850] = 0e0;
		v[851] = 0e0;
		b852 = b266;
		if (b852) {
			v[848] = v[400] * v[814] + v[848];
			v[851] = v[716] * v[814];
			v[848] = v[2464] * v[841] + v[848];
			v[850] = v[271] * v[716] * v[841];
		}
		else {
		};
		v[853] = 0e0;
		v[854] = 0e0;
		b855 = b266;
		if (b855) {
			v[853] = -(v[400] * v[814]);
			v[851] = -(v[713] * v[814]) + v[851];
			v[853] = v[2465] * v[841] + v[853];
			v[854] = v[271] * v[713] * v[841];
			v[841] = 0e0;
		}
		else {
		};
		b856 = b273;
		if (b856) {
			v[853] = v[568] * v[814] + v[853];
			v[814] = 0e0;
		}
		else {
		};
		v[2530] = v[2527] + v[325] * v[821] + v[328] * v[827] + v[331] * v[832] - v[334] * v[842] - v[337] * v[848]
			- v[340] * v[853];
		v[973] = i2507 * v[200] + i2510 * v[205] + i2513 * v[211] - i2516 * v[302] - i2519 * v[307] - i2522 * v[313] + v[119] * v[821]
			+ v[122] * v[827] + v[125] * v[832] - v[221] * v[842] - v[224] * v[848] - v[227] * v[853];
		v[974] = i2508 * v[200] + i2511 * v[205] + i2514 * v[211] - i2517 * v[302] - i2520 * v[307] - i2523 * v[313] + v[120] * v[821]
			+ v[123] * v[827] + v[126] * v[832] - v[222] * v[842] - v[225] * v[848] - v[228] * v[853];
		v[975] = i2509 * v[200] + i2512 * v[205] + i2515 * v[211] - i2518 * v[302] - i2521 * v[307] - i2524 * v[313] + v[121] * v[821]
			+ v[124] * v[827] + v[127] * v[832] - v[223] * v[842] - v[226] * v[848] - v[229] * v[853];
		v[2546] = v[341] * v[973] + v[342] * v[974] + v[343] * v[975];
		v[982] = -(v[324] * v[821]) - v[327] * v[827] - v[330] * v[832] + v[333] * v[842] + v[336] * v[848] + v[339] * v[853]
			+ v[872] + v[873] + v[874] - v[875] - v[876] - v[877];
		v[984] = -(v[323] * v[821]) - v[326] * v[827] - v[329] * v[832] + v[332] * v[842] + v[335] * v[848] + v[338] * v[853]
			+ v[880] + v[881] + v[882] - v[883] - v[884] - v[885];
		v[992] = v[2546] / v[351];
		v[2525] = v[356] * v[992];
		v[993] = v[2525] * v[341] + v[357] * v[973];
		v[2526] = -(v[359] * v[993]);
		v[994] = v[2525] * v[342] + v[357] * v[974];
		v[995] = v[2525] * v[343] + v[357] * v[975];
		v[996] = v[358] * v[992] + v[350] * v[993];
		v[2541] = v[13] * v[996];
		v[1060] = -(v[2541] * v[520]);
		v[997] = v[2526] * v[313] - v[358] * (v[359] * v[853] + v[313] * v[994]);
		v[998] = v[2526] * v[307] - v[358] * (v[359] * v[848] + v[307] * v[994]);
		v[999] = v[2526] * v[302] - v[358] * (v[359] * v[842] + v[302] * v[994]);
		v[1000] = -(v[211] * v[2526]) + v[358] * (v[359] * v[832] + v[211] * v[994]);
		v[1001] = -(v[205] * v[2526]) + v[358] * (v[359] * v[827] + v[205] * v[994]);
		v[1002] = -(v[200] * v[2526]) + v[358] * (v[359] * v[821] + v[200] * v[994]);
		v[1003] = v[359] * v[992] + v[350] * v[994];
		v[2538] = v[1003] * v[13];
		v[2529] = v[1003] * v[359] + v[358] * v[996];
		v[1061] = -(v[2538] * v[528]);
		v[1004] = v[358] * v[993] + v[359] * v[994];
		v[2543] = (v[2527] * v[536] + v[821] * v[869] + v[827] * v[900] + v[832] * v[922] + v[842] * v[941] + v[848] * v[960]
			+ v[853] * v[979] + v[1005] * v[993] + v[1006] * v[994]) / 2e0;
		v[1009] = v[2530] * v[360] + v[1008] * v[995];
		v[1032] = v[1009] * v[2472];
		v[1010] = v[360] * v[992] + v[350] * v[995];
		v[2528] = v[1010] * v[13];
		v[1062] = -(v[2528] * v[536]);
		v[1011] = v[1002] * v[323] + v[1001] * v[326] + v[1000] * v[329] + v[338] * v[997] + v[335] * v[998] + v[332] * v[999];
		v[1012] = v[1002] * v[324] + v[1001] * v[327] + v[1000] * v[330] + v[339] * v[997] + v[336] * v[998] + v[333] * v[999];
		v[1013] = v[2528] * v[360] + v[705] * v[995];
		v[2570] = v[1061] - v[1013] * v[359];
		v[2569] = v[1060] - v[1013] * v[358];
		v[1015] = v[15] * (i2508 * v[751] + i2507 * v[752]) + v[13] * (v[1003] * v[323] + v[324] * v[996]);
		v[1016] = v[15] * (i2511 * v[751] + i2510 * v[752]) + v[13] * (v[1003] * v[326] + v[327] * v[996]);
		v[1017] = v[15] * (i2514 * v[751] + i2513 * v[752]) + v[13] * (v[1003] * v[329] + v[330] * v[996]);
		v[1018] = v[15] * (i2517 * v[751] + i2516 * v[752]) + v[13] * (v[1003] * v[332] + v[333] * v[996]);
		v[1019] = v[15] * (i2520 * v[751] + i2519 * v[752]) + v[13] * (v[1003] * v[335] + v[336] * v[996]);
		v[1020] = v[15] * (i2523 * v[751] + i2522 * v[752]) + v[13] * (v[1003] * v[338] + v[339] * v[996]);
		v[2540] = v[1015] * v[200] + v[1016] * v[205] + v[1017] * v[211] - v[1018] * v[302] - v[1019] * v[307] - v[1020] * v[313];
		v[1021] = v[1004] * v[2472] + v[13] * v[2529];
		v[2571] = -(v[1021] * v[360]);
		v[2572] = -v[1062] - v[2571] + v[2697] * v[995];
		v[2556] = v[1062] + v[2571];
		v[1025] = v[1008] * v[1021] + v[1010] * v[1920] + v[2527] * v[669] + v[2703] * v[821] + v[2702] * v[827]
			+ v[2701] * v[832] + v[2700] * v[842] + v[2699] * v[848] + v[2698] * v[853] + v[587] * v[992] + v[705] * (v[1022] * v[993]
				+ v[1023] * v[994]) + v[703] * (v[2543] + (v[1024] * v[995]) / 2e0) + v[2582] * (v[2530] * v[705] + v[2539] * v[751]
					+ v[2537] * v[752] + v[1008] * (v[2528] + v[703] * v[995]));
		v[2544] = (v[15] * (i2507 * v[522] + i2510 * v[523] + i2513 * v[524] + i2516 * v[525] + i2519 * v[526] + i2522 * v[527])
			+ v[2537] * v[536] - v[528] * v[982]) / 2e0;
		v[1028] = v[1013] * v[1023] + v[1003] * v[1031] + v[2540] * v[358] + (v[1011] / 2e0 + v[2544]) * v[703] + v[2567] * v[821]
			+ v[2564] * v[827] + v[2561] * v[832] + v[2558] * v[842] + v[2553] * v[848] + v[2549] * v[853] + v[586] * v[992]
			+ v[1030] * v[993] + v[706] * v[994] + v[2583] * (v[1032] - v[752] * v[982] + v[1023] * (v[2538] + v[703] * v[994]));
		v[2545] = (v[15] * (i2508 * v[522] + i2511 * v[523] + i2514 * v[524] + i2517 * v[525] + i2520 * v[526] + i2523 * v[527])
			+ v[2539] * v[536] - v[520] * v[984]) / 2e0;
		v[1033] = v[1013] * v[1022] + v[2540] * v[359] + (v[1012] / 2e0 + v[2545]) * v[703] + v[2566] * v[821] + v[2563] * v[827]
			+ v[2560] * v[832] + v[2557] * v[842] + v[2551] * v[848] + v[2548] * v[853] + v[585] * v[992] + v[706] * v[993] + v[2584] *
			(v[1032] - v[751] * v[984] + v[1022] * (v[2541] + v[703] * v[993])) + v[1030] * v[994] + v[1031] * v[996];
		v[1039] = v[2649] * (v[2529] + v[1010] * v[360]) - v[1034] * v[2546] * v[708] + v[356] * (v[1033] * v[341]
			+ v[1028] * v[342] + v[1025] * v[343] + v[707] * v[973] + v[704] * v[974] + v[702] * v[975]) + v[1035] * v[355] * v[992]
			+ v[585] * v[993] + v[586] * v[994] + v[2547] * v[995] + v[13] * (v[1009] * v[2542] + v[1012] * v[358] + v[1011] * v[359]
				+ 2e0 * ((v[1004] * v[1027]) / 2e0 + v[2545] * v[358] + v[2544] * v[359] + v[2543] * v[360] + v[2499] * v[993]
					+ v[2498] * v[994] + v[2497] * v[995]));
		v[1041] = v[1025] * v[357] + v[2525] * v[702] + (v[1039] * v[343] + v[708] * v[975]) / v[351];
		v[1043] = v[1028] * v[357] + v[2525] * v[704] + (v[1039] * v[342] + v[708] * v[974]) / v[351];
		v[1045] = v[1033] * v[357] + v[2525] * v[707] + (v[1039] * v[341] + v[708] * v[973]) / v[351];
		v[1051] = i2524 * v[1052] + i2523 * v[1053] + i2522 * v[1054] - v[1045] * v[227] - v[1043] * v[228] - v[1041] * v[229]
			+ v[1020] * v[2479] + v[1013] * v[2531] + v[1060] * v[338] + v[1061] * v[339] + v[2556] * v[340] + v[2705] * v[993]
			+ v[2706] * v[994] + v[2704] * v[995];
		v[1063] = i2521 * v[1052] + i2520 * v[1053] + i2519 * v[1054] - v[1045] * v[224] - v[1043] * v[225] - v[1041] * v[226]
			+ v[1019] * v[2479] + v[1013] * v[2532] + v[1060] * v[335] + v[1061] * v[336] + v[2556] * v[337] + v[2707] * v[993]
			+ v[2708] * v[994] + v[2709] * v[995];
		v[1069] = i2518 * v[1052] + i2517 * v[1053] + i2516 * v[1054] - v[1045] * v[221] - v[1043] * v[222] - v[1041] * v[223]
			+ v[1018] * v[2479] + v[1013] * v[2533] + v[1060] * v[332] + v[1061] * v[333] + v[2556] * v[334] + v[2710] * v[993]
			+ v[2711] * v[994] + v[2712] * v[995];
		v[1075] = -(i2515 * v[1052]) - i2514 * v[1053] - i2513 * v[1054] + v[1045] * v[125] + v[1043] * v[126] + v[1041] * v[127]
			- v[1017] * v[2479] + v[1013] * v[2534] - v[1060] * v[329] - v[1061] * v[330] - v[2556] * v[331] + v[2713] * v[993]
			+ v[2714] * v[994] + v[2715] * v[995];
		v[1081] = -(i2512 * v[1052]) - i2511 * v[1053] - i2510 * v[1054] + v[1045] * v[122] + v[1043] * v[123] + v[1041] * v[124]
			- v[1016] * v[2479] + v[1013] * v[2535] - v[1060] * v[326] - v[1061] * v[327] - v[2556] * v[328] + v[2716] * v[993]
			+ v[2717] * v[994] + v[2718] * v[995];
		v[1087] = -(i2509 * v[1052]) - i2508 * v[1053] - i2507 * v[1054] + v[1045] * v[119] + v[1043] * v[120] + v[1041] * v[121]
			- v[1015] * v[2479] + v[1013] * v[2536] - v[1060] * v[323] - v[1061] * v[324] - v[2556] * v[325] + v[2719] * v[993]
			+ v[2720] * v[994] + v[2721] * v[995];
		v[1088] = 0e0;
		b1089 = b273;
		if (b1089) {
			v[1088] = v[1051] * v[568];
		}
		else {
		};
		v[1090] = 0e0;
		b1091 = b266;
		if (b1091) {
			v[854] = -(v[1051] * v[400]) + v[854];
			v[1090] = v[1051] * v[2465];
			v[1088] = v[1088] + v[854];
		}
		else {
		};
		b1092 = b266;
		if (b1092) {
			v[850] = -(v[1063] * v[400]) + v[850];
			v[1088] = v[1088] - v[850];
			v[1090] = v[1090] + v[271] * (v[1063] * v[397] - v[851]);
		}
		else {
		};
		v[1093] = 0e0;
		b1094 = b259;
		if (b1094) {
			v[847] = -(v[1063] * v[394]) + v[847];
			v[1093] = v[1063] * v[2463];
			v[1088] = v[1088] + v[847];
		}
		else {
		};
		b1095 = b259;
		if (b1095) {
			v[845] = -(v[1069] * v[394]) + v[845];
			v[1093] = v[1093] + v[264] * (v[1069] * v[730] - v[844]);
			v[1088] = v[1088] - v[845];
		}
		else {
		};
		b1096 = b256;
		if (b1096) {
			v[1088] = v[1088] + v[1069] * v[558];
		}
		else {
		};
		v[1097] = v[1090];
		b1098 = b239;
		if (b1098) {
			v[1088] = v[1088] + v[1097] * v[244];
		}
		else {
		};
		v[1099] = v[1093];
		b1100 = b239;
		if (b1100) {
			v[1088] = v[1088] - v[1099] * v[244];
		}
		else {
		};
		v[1101] = 0e0;
		b1102 = b171;
		if (b1102) {
			v[1101] = v[1075] * v[550];
		}
		else {
		};
		v[1103] = 0e0;
		b1104 = b164;
		if (b1104) {
			v[833] = -(v[1075] * v[381]) + v[833];
			v[1103] = v[1075] * v[2459];
			v[1101] = v[1101] + v[833];
		}
		else {
		};
		b1105 = b164;
		if (b1105) {
			v[829] = -(v[1081] * v[381]) + v[829];
			v[1101] = v[1101] - v[829];
			v[1103] = v[1103] + v[169] * (v[1081] * v[378] - v[830]);
		}
		else {
		};
		v[1106] = 0e0;
		b1107 = b157;
		if (b1107) {
			v[826] = -(v[1081] * v[375]) + v[826];
			v[1106] = v[1081] * v[2457];
			v[1101] = v[1101] + v[826];
		}
		else {
		};
		b1108 = b157;
		if (b1108) {
			v[824] = -(v[1087] * v[375]) + v[824];
			v[1106] = v[1106] + v[162] * (v[1087] * v[745] - v[823]);
			v[1101] = v[1101] - v[824];
		}
		else {
		};
		b1109 = b154;
		if (b1109) {
			v[1101] = v[1101] + v[1087] * v[540];
		}
		else {
		};
		v[1110] = v[1103];
		b1111 = b137;
		if (b1111) {
			v[1101] = v[1101] + v[1110] * v[142];
		}
		else {
		};
		v[1112] = v[1106];
		b1113 = b137;
		if (b1113) {
			v[1101] = v[1101] - v[1112] * v[142];
		}
		else {
		};
		v[4402] = v[1045] * v[200] + v[1101] * v[444] + v[1088] * v[462] + v[712] * v[821] + v[15] * (-(v[200] * v[2569])
			+ v[2538] * v[522] + v[1002] * v[752] + v[821] * v[911] + v[2728] * v[993]);
		v[4403] = v[1043] * v[200] + v[1101] * v[445] + v[1088] * v[463] + v[711] * v[821] + v[15] * (-(v[200] * v[2570])
			+ v[2541] * v[522] + v[1002] * v[751] + v[821] * v[909] + v[2729] * v[994]);
		v[4404] = v[1041] * v[200] + v[1101] * v[446] + v[1088] * v[464] + v[709] * v[821] + v[15] * (v[200] * v[2572]
			+ v[821] * v[907]);
		v[4405] = v[1045] * v[205] + v[1101] * v[447] + v[1088] * v[465] + v[712] * v[827] + v[15] * (-(v[205] * v[2569])
			+ v[2538] * v[523] + v[1001] * v[752] + v[827] * v[911] + v[2730] * v[993]);
		v[4406] = v[1043] * v[205] + v[1101] * v[448] + v[1088] * v[466] + v[711] * v[827] + v[15] * (-(v[205] * v[2570])
			+ v[2541] * v[523] + v[1001] * v[751] + v[827] * v[909] + v[2731] * v[994]);
		v[4407] = v[1041] * v[205] + v[1101] * v[449] + v[1088] * v[467] + v[709] * v[827] + v[15] * (v[205] * v[2572]
			+ v[827] * v[907]);
		v[4408] = v[1045] * v[211] + v[1101] * v[450] + v[1088] * v[468] + v[712] * v[832] + v[15] * (-(v[211] * v[2569])
			+ v[2538] * v[524] + v[1000] * v[752] + v[832] * v[911] + v[2732] * v[993]);
		v[4409] = v[1043] * v[211] + v[1101] * v[451] + v[1088] * v[469] + v[711] * v[832] + v[15] * (-(v[211] * v[2570])
			+ v[2541] * v[524] + v[1000] * v[751] + v[832] * v[909] + v[2733] * v[994]);
		v[4410] = v[1041] * v[211] + v[1101] * v[452] + v[1088] * v[470] + v[709] * v[832] + v[15] * (v[211] * v[2572]
			+ v[832] * v[907]);
		v[4411] = -(v[1045] * v[302]) + v[1101] * v[453] + v[1088] * v[471] - v[712] * v[842] + v[15] * (v[2569] * v[302]
			+ v[2538] * v[525] - v[842] * v[911] + v[2734] * v[993] + v[752] * v[999]);
		v[4412] = -(v[1043] * v[302]) + v[1101] * v[454] + v[1088] * v[472] - v[711] * v[842] + v[15] * (v[2570] * v[302]
			+ v[2541] * v[525] - v[842] * v[909] + v[2735] * v[994] + v[751] * v[999]);
		v[4413] = -(v[1041] * v[302]) + v[1101] * v[455] + v[1088] * v[473] - v[709] * v[842] + v[15] * (-(v[2572] * v[302])
			- v[842] * v[907]);
		v[4414] = -(v[1045] * v[307]) + v[1101] * v[456] + v[1088] * v[474] - v[712] * v[848] + v[15] * (v[2569] * v[307]
			+ v[2538] * v[526] - v[848] * v[911] + v[2736] * v[993] + v[752] * v[998]);
		v[4415] = -(v[1043] * v[307]) + v[1101] * v[457] + v[1088] * v[475] - v[711] * v[848] + v[15] * (v[2570] * v[307]
			+ v[2541] * v[526] - v[848] * v[909] + v[2737] * v[994] + v[751] * v[998]);
		v[4416] = -(v[1041] * v[307]) + v[1101] * v[458] + v[1088] * v[476] - v[709] * v[848] + v[15] * (-(v[2572] * v[307])
			- v[848] * v[907]);
		v[4417] = -(v[1045] * v[313]) + v[1101] * v[459] + v[1088] * v[477] - v[712] * v[853] + v[15] * (v[2569] * v[313]
			+ v[2538] * v[527] - v[853] * v[911] + v[2738] * v[993] + v[752] * v[997]);
		v[4418] = -(v[1043] * v[313]) + v[1101] * v[460] + v[1088] * v[478] - v[711] * v[853] + v[15] * (v[2570] * v[313]
			+ v[2541] * v[527] - v[853] * v[909] + v[2739] * v[994] + v[751] * v[997]);
		v[4419] = -(v[1041] * v[313]) + v[1101] * v[461] + v[1088] * v[479] - v[709] * v[853] + v[15] * (-(v[2572] * v[313])
			- v[853] * v[907]);
		Rc[i667 - 1] += v[4055 + i667];
		for (i792 = i667; i792 <= 18; i792++) {
			v[1157] = v[4401 + i792];
			Kc[i667 - 1][i792 - 1] += v[1157];
			if (i667 != i792) {
				Kc[i792 - 1][i667 - 1] += v[1157];
			}
			else {
			};
		};/* end for */
	};/* end for */
	v[1215] = 0e0;
	v[1216] = 0e0;
	v[1217] = 0e0;
	b1218 = b22;
	if (b1218) {
		b1219 = b598;
		if (b1219) {
			v[1217] = 0e0;
			v[1216] = 0e0;
			v[1215] = 0e0;
		}
		else {
		};
	}
	else {
		b1220 = b632;
		if (b1220) {
			v[1217] = 0e0;
			v[1216] = 0e0;
			v[1215] = 0e0;
		}
		else {
		};
	};
	v[2575] = v[1215] * v[14];
	v[2574] = v[1216] * v[14];
	v[2573] = v[1217] * v[14];
	v[1261] = v[1300] * v[2573];
	v[1262] = v[1302] * v[2573];
	v[1301] = v[1300] * v[2574];
	v[1303] = v[1302] * v[2574];
	v[1341] = v[1300] * v[2575];
	v[1342] = v[1302] * v[2575];
	v[1343] = v[10] * v[1217];
	v[1696] = -(v[1343] * v[360]);
	v[1344] = v[10] * v[1216];
	v[1700] = -(v[1344] * v[359]);
	v[1702] = v[1696] + v[1700];
	v[1345] = v[10] * v[1215];
	v[1701] = -(v[1345] * v[358]);
	v[1704] = v[1700] + v[1701];
	v[1703] = v[1696] + v[1701];
	v[1346] = 0e0;
	v[1347] = 0e0;
	v[1348] = 0e0;
	v[1349] = 0e0;
	v[1350] = 0e0;
	v[1351] = 0e0;
	v[1352] = 0e0;
	v[1353] = 0e0;
	v[1354] = 0e0;
	v[1355] = 0e0;
	v[1356] = 0e0;
	v[1357] = 0e0;
	b1358 = b23;
	if (b1358) {
		v[1359] = -(v[1343] * v[517]);
		v[1360] = -(v[1344] * v[516]);
		v[1362] = v[1345] * v[1361] + v[1702] * v[358];
		v[1364] = v[1344] * v[1363] + v[1703] * v[359];
		v[1366] = v[1343] * v[1365] + v[1704] * v[360];
		v[1348] = v[1343] * v[2576] + v[1704] * v[517];
		v[1347] = v[1344] * v[2577] + v[1703] * v[516];
		v[1368] = -(v[1345] * v[515]);
		v[1346] = -(v[1345] * v[2578]) + v[1702] * v[515];
		v[1349] = v[1362] * v[5];
		v[1350] = v[1362] * v[6];
		v[1351] = v[1362] * v[7];
		v[1371] = -(v[1362] * v[276]);
		v[1372] = -(v[1362] * v[269]);
		v[1373] = -(v[1362] * v[262]);
		v[1374] = v[1362] * v[174];
		v[1375] = v[1362] * v[167];
		v[1376] = v[1362] * v[160];
		v[1352] = v[1364] * v[5];
		v[1353] = v[1364] * v[6];
		v[1354] = v[1364] * v[7];
		v[1377] = -(v[1364] * v[276]);
		v[1378] = -(v[1364] * v[269]);
		v[1379] = -(v[1364] * v[262]);
		v[1380] = v[1364] * v[174];
		v[1381] = v[1364] * v[167];
		v[1382] = v[1364] * v[160];
		v[1355] = v[1366] * v[5];
		v[1356] = v[1366] * v[6];
		v[1357] = v[1366] * v[7];
		v[1383] = -(v[1366] * v[276]);
		v[1384] = -(v[1366] * v[269]);
		v[1385] = -(v[1366] * v[262]);
		v[1386] = v[1366] * v[174];
		v[1387] = v[1366] * v[167];
		v[1388] = v[1366] * v[160];
	}
	else {
		v[1376] = 0e0;
		v[1382] = 0e0;
		v[1388] = 0e0;
		v[1375] = 0e0;
		v[1381] = 0e0;
		v[1387] = 0e0;
		v[1374] = 0e0;
		v[1380] = 0e0;
		v[1386] = 0e0;
		v[1373] = 0e0;
		v[1379] = 0e0;
		v[1385] = 0e0;
		v[1372] = 0e0;
		v[1378] = 0e0;
		v[1384] = 0e0;
		v[1371] = 0e0;
		v[1377] = 0e0;
		v[1383] = 0e0;
		v[1368] = 0e0;
		v[1360] = 0e0;
		v[1359] = 0e0;
	};
	v[2611] = v[1349] / 2e0;
	v[2612] = v[1353] / 2e0;
	v[2613] = v[1357] / 2e0;
	b1389 = b480;
	if (b1389) {
		v[1423] = -(v[1350] * v[497]);
		v[1418] = v[1351] * v[497];
		v[1405] = v[1354] * v[497];
		v[1392] = -(v[1357] * v[2489]);
		v[1396] = v[1356] * v[497];
		v[1400] = v[1355] * v[497];
		v[1404] = v[1396] + v[1405];
		v[1409] = -(v[1353] * v[2489]);
		v[1413] = v[1352] * v[497];
		v[1417] = v[1400] + v[1418];
		v[1424] = v[1413] - v[1423];
		v[1426] = v[1356] * v[1394] + v[1355] * v[1398] + v[1354] * v[1402] + v[1352] * v[1411] + v[1351] * v[1415]
			+ v[1350] * v[1420] + v[1425] * v[2611] + v[1407] * v[2612] + v[1390] * v[2613];
		v[1650] = v[1392] + v[1409] - (4e0 * v[1426]) / (v[1430] * v[1430]);
		v[2610] = 4e0 * v[1650];
		v[1648] = -v[1392] + v[1650] - v[1349] * v[2489];
		v[2609] = 4e0 * (v[1392] - v[1409] + v[1648]);
		v[1431] = v[1413] + v[1423] + v[1404] * v[2488] + v[1417] * v[2579] + 2e0 * v[1648] * v[496];
		v[1433] = (-2e0 * v[1400] + 2e0 * v[1418] + v[1424] * v[494] + v[2609] * v[495] + v[1404] * v[496]) / 2e0;
		v[1434] = (2e0 * v[1396] - 2e0 * v[1405] + v[2610] * v[494] + v[1424] * v[495] + v[1417] * v[496]) / 2e0;
		v[2580] = v[1434] * v[482] + v[1433] * v[483] + v[1431] * v[484];
		v[1636] = v[2580] * v[493];
		v[1633] = v[2580] * v[487];
		v[1437] = v[1633] * v[492] + v[1636] / (Power(cos(v[1435]), 2) * sqrt(v[1637]));
		v[2581] = v[1437] / v[485];
		v[1438] = v[1431] * v[2487] + v[2581] * v[484];
		v[1440] = v[1433] * v[2487] + v[2581] * v[483];
		v[1441] = v[1434] * v[2487] + v[2581] * v[482];
		v[1346] = v[1346] - v[1438] * v[368] + v[1440] * v[369];
		v[1348] = v[1348] - v[1440] * v[367] + v[1441] * v[368];
		v[1347] = v[1347] + v[1438] * v[367] - v[1441] * v[369];
	}
	else {
	};
	v[1348] = v[1348] + v[1359] * v[2582];
	v[1347] = v[1347] + v[1360] * v[2583];
	v[1346] = v[1346] + v[1368] * v[2584];
	v[1446] = v[1346] * v[341] + v[1347] * v[342] + v[1348] * v[343];
	v[1448] = v[1446] * v[356];
	v[2585] = v[1448] / v[351];
	v[1449] = v[2585] * v[343] + v[1348] * v[357];
	v[1450] = v[2585] * v[342] + v[1347] * v[357];
	v[1451] = v[2585] * v[341] + v[1346] * v[357];
	v[1371] = v[1371] - v[1451] * v[313] + v[1341] * v[570];
	v[1372] = v[1372] - v[1451] * v[307] + v[1341] * v[565];
	v[1373] = v[1373] - v[1451] * v[302] + v[1341] * v[562];
	v[1377] = v[1377] - v[1450] * v[313] + v[1301] * v[570];
	v[1378] = v[1378] - v[1450] * v[307] + v[1301] * v[565];
	v[1379] = v[1379] - v[1450] * v[302] + v[1301] * v[562];
	v[1383] = v[1383] - v[1449] * v[313] + v[1261] * v[570];
	v[1384] = v[1384] - v[1449] * v[307] + v[1261] * v[565];
	v[1385] = v[1385] - v[1449] * v[302] + v[1261] * v[562];
	b1460 = b273;
	b1463 = b266;
	b1464 = b266;
	b1467 = b259;
	b1469 = b259;
	b1470 = b256;
	b1472 = b239;
	b1474 = b239;
	v[1374] = v[1374] + v[1451] * v[211] + v[1342] * v[552];
	v[1375] = v[1375] + v[1451] * v[205] + v[1342] * v[547];
	v[1376] = v[1376] + v[1451] * v[200] + v[1342] * v[544];
	v[1380] = v[1380] + v[1450] * v[211] + v[1303] * v[552];
	v[1381] = v[1381] + v[1450] * v[205] + v[1303] * v[547];
	v[1382] = v[1382] + v[1450] * v[200] + v[1303] * v[544];
	v[1386] = v[1386] + v[1449] * v[211] + v[1262] * v[552];
	v[1387] = v[1387] + v[1449] * v[205] + v[1262] * v[547];
	v[1388] = v[1388] + v[1449] * v[200] + v[1262] * v[544];
	b1486 = b171;
	b1489 = b164;
	b1490 = b164;
	b1493 = b157;
	b1495 = b157;
	b1496 = b154;
	b1498 = b137;
	b1500 = b137;
	v[4424] = v[1376] + (-(v[1217] * v[1259]) - v[1216] * v[1298] - v[1215] * v[1339]) * v[2586] - v[1707] * v[600]
		- v[1726] * v[601] - v[1745] * v[602];
	v[4425] = v[1382] + (-(v[1217] * v[1257]) - v[1216] * v[1296] - v[1215] * v[1337]) * v[2586] - v[1708] * v[600]
		- v[1727] * v[601] - v[1746] * v[602];
	v[4426] = v[1388] + (-(v[1217] * v[1255]) - v[1216] * v[1294] - v[1215] * v[1335]) * v[2586] - v[1709] * v[600]
		- v[1728] * v[601] - v[1747] * v[602];
	v[4427] = v[1375] + (-(v[1217] * v[1253]) - v[1216] * v[1292] - v[1215] * v[1333]) * v[2586] - v[1710] * v[600]
		- v[1729] * v[601] - v[1748] * v[602];
	v[4428] = v[1381] + (-(v[1217] * v[1251]) - v[1216] * v[1290] - v[1215] * v[1331]) * v[2586] - v[1711] * v[600]
		- v[1730] * v[601] - v[1749] * v[602];
	v[4429] = v[1387] + (-(v[1217] * v[1249]) - v[1216] * v[1288] - v[1215] * v[1329]) * v[2586] - v[1712] * v[600]
		- v[1731] * v[601] - v[1750] * v[602];
	v[4430] = v[1374] + (-(v[1217] * v[1247]) - v[1216] * v[1286] - v[1215] * v[1327]) * v[2586] - v[1713] * v[600]
		- v[1732] * v[601] - v[1751] * v[602];
	v[4431] = v[1380] + (-(v[1217] * v[1245]) - v[1216] * v[1284] - v[1215] * v[1325]) * v[2586] - v[1714] * v[600]
		- v[1733] * v[601] - v[1752] * v[602];
	v[4432] = v[1386] + (-(v[1217] * v[1243]) - v[1216] * v[1282] - v[1215] * v[1323]) * v[2586] - v[1715] * v[600]
		- v[1734] * v[601] - v[1753] * v[602];
	v[4433] = v[1373] + (-(v[1217] * v[1241]) - v[1216] * v[1280] - v[1215] * v[1321]) * v[2586] - v[1716] * v[600]
		- v[1735] * v[601] - v[1754] * v[602];
	v[4434] = v[1379] + (-(v[1217] * v[1239]) - v[1216] * v[1278] - v[1215] * v[1319]) * v[2586] - v[1717] * v[600]
		- v[1736] * v[601] - v[1755] * v[602];
	v[4435] = v[1385] + (-(v[1217] * v[1237]) - v[1216] * v[1276] - v[1215] * v[1317]) * v[2586] - v[1718] * v[600]
		- v[1737] * v[601] - v[1756] * v[602];
	v[4436] = v[1372] + (-(v[1217] * v[1235]) - v[1216] * v[1274] - v[1215] * v[1315]) * v[2586] - v[1719] * v[600]
		- v[1738] * v[601] - v[1757] * v[602];
	v[4437] = v[1378] + (-(v[1217] * v[1233]) - v[1216] * v[1272] - v[1215] * v[1313]) * v[2586] - v[1720] * v[600]
		- v[1739] * v[601] - v[1758] * v[602];
	v[4438] = v[1384] + (-(v[1217] * v[1231]) - v[1216] * v[1270] - v[1215] * v[1311]) * v[2586] - v[1721] * v[600]
		- v[1740] * v[601] - v[1759] * v[602];
	v[4439] = v[1371] + (-(v[1217] * v[1229]) - v[1216] * v[1268] - v[1215] * v[1309]) * v[2586] - v[1722] * v[600]
		- v[1741] * v[601] - v[1760] * v[602];
	v[4440] = v[1377] + (-(v[1217] * v[1227]) - v[1216] * v[1266] - v[1215] * v[1307]) * v[2586] - v[1723] * v[600]
		- v[1742] * v[601] - v[1761] * v[602];
	v[4441] = v[1383] + (-(v[1217] * v[1225]) - v[1216] * v[1264] - v[1215] * v[1305]) * v[2586] - v[1724] * v[600]
		- v[1743] * v[601] - v[1762] * v[602];
	for (i1213 = 1; i1213 <= 18; i1213++) {
		i2604 = (i1213 == 16 ? 1 : 0);
		i2603 = (i1213 == 13 ? 1 : 0);
		i2602 = (i1213 == 10 ? 1 : 0);
		i2601 = (i1213 == 7 ? 1 : 0);
		i2600 = (i1213 == 4 ? 1 : 0);
		i2599 = (i1213 == 1 ? 1 : 0);
		i2598 = (i1213 == 17 ? 1 : 0);
		i2597 = (i1213 == 14 ? 1 : 0);
		i2596 = (i1213 == 11 ? 1 : 0);
		i2595 = (i1213 == 8 ? 1 : 0);
		i2594 = (i1213 == 5 ? 1 : 0);
		i2593 = (i1213 == 2 ? 1 : 0);
		i2592 = (i1213 == 18 ? 1 : 0);
		i2591 = (i1213 == 15 ? 1 : 0);
		i2590 = (i1213 == 12 ? 1 : 0);
		i2589 = (i1213 == 9 ? 1 : 0);
		i2588 = (i1213 == 6 ? 1 : 0);
		i2587 = (i1213 == 3 ? 1 : 0);
		v[2624] = i2599 * v[444] + i2593 * v[445] + i2587 * v[446] + i2600 * v[447] + i2594 * v[448] + i2588 * v[449] + i2601 * v[450]
			+ i2595 * v[451] + i2589 * v[452] + i2602 * v[453] + i2596 * v[454] + i2590 * v[455] + i2603 * v[456] + i2597 * v[457]
			+ i2591 * v[458] + i2604 * v[459] + i2598 * v[460] + i2592 * v[461];
		v[2623] = i2599 * v[462] + i2593 * v[463] + i2587 * v[464] + i2600 * v[465] + i2594 * v[466] + i2588 * v[467] + i2601 * v[468]
			+ i2595 * v[469] + i2589 * v[470] + i2602 * v[471] + i2596 * v[472] + i2590 * v[473] + i2603 * v[474] + i2597 * v[475]
			+ i2591 * v[476] + i2604 * v[477] + i2598 * v[478] + i2592 * v[479];
		v[1533] = i2599 * v[15];
		v[1534] = i2593 * v[15];
		v[1535] = i2587 * v[15];
		v[1536] = i2600 * v[15];
		v[1537] = i2594 * v[15];
		v[1538] = i2588 * v[15];
		v[1539] = i2601 * v[15];
		v[1540] = i2595 * v[15];
		v[1541] = i2589 * v[15];
		v[1551] = i2602 * v[15];
		v[1552] = i2596 * v[15];
		v[1553] = i2590 * v[15];
		v[1554] = i2603 * v[15];
		v[1555] = i2597 * v[15];
		v[1556] = i2591 * v[15];
		v[1557] = i2604 * v[15];
		v[1558] = i2598 * v[15];
		v[1559] = i2592 * v[15];
		v[1562] = i2587 * v[200] + i2588 * v[205] + i2589 * v[211] - i2590 * v[302] - i2591 * v[307] - i2592 * v[313];
		v[1566] = i2593 * v[200] + i2594 * v[205] + i2595 * v[211] - i2596 * v[302] - i2597 * v[307] - i2598 * v[313];
		v[1570] = i2599 * v[200] + i2600 * v[205] + i2601 * v[211] - i2602 * v[302] - i2603 * v[307] - i2604 * v[313];
		v[2648] = v[1570] * v[341] + v[1566] * v[342] + v[1562] * v[343];
		v[1572] = v[2648] / v[351];
		v[1573] = v[1572] * v[356];
		v[1583] = v[1573] * v[343] + v[1562] * v[357];
		v[2620] = 2e0 * v[1583];
		v[1579] = v[1573] * v[342] + v[1566] * v[357];
		v[2618] = 2e0 * v[1579];
		v[1575] = v[1573] * v[341] + v[1570] * v[357];
		v[2617] = 2e0 * v[1575];
		v[1576] = v[1368] * v[2617];
		v[1577] = v[1575];
		v[1580] = v[1360] * v[2618];
		v[1581] = v[1579];
		v[1584] = v[1359] * v[2620];
		v[1585] = v[1583];
		v[1586] = 0e0;
		v[1587] = 0e0;
		v[1588] = 0e0;
		v[1589] = 0e0;
		v[1590] = 0e0;
		v[1591] = 0e0;
		v[1592] = 0e0;
		v[1593] = 0e0;
		v[1594] = 0e0;
		v[1595] = 0e0;
		v[1596] = 0e0;
		v[1597] = 0e0;
		v[1598] = 0e0;
		v[1599] = 0e0;
		v[1600] = 0e0;
		v[1601] = 0e0;
		v[1602] = 0e0;
		v[1603] = 0e0;
		v[1604] = 0e0;
		v[1605] = 0e0;
		v[1606] = 0e0;
		v[1607] = 0e0;
		v[1608] = 0e0;
		v[1609] = 0e0;
		v[1610] = 0e0;
		v[1611] = 0e0;
		v[1612] = 0e0;
		v[1613] = 0e0;
		v[1614] = 0e0;
		v[1615] = 0e0;
		v[1616] = 0e0;
		v[1617] = 0e0;
		b1618 = b480;
		if (b1618) {
			v[1619] = -(v[1581] * v[369]);
			v[1620] = v[1581] * v[367];
			v[1621] = v[1619] + v[1585] * v[368];
			v[1622] = -(v[1585] * v[367]);
			v[1623] = v[1622] + v[1577] * v[369];
			v[1624] = v[1620] - v[1577] * v[368];
			v[2608] = v[1434] * v[1621] + v[1433] * v[1623] + v[1431] * v[1624];
			v[2606] = v[1621] * v[482] + v[1623] * v[483] + v[1624] * v[484];
			v[1625] = v[2606] / v[485];
			v[2607] = v[1625] * v[2772];
			v[1635] = v[1625] * v[2605];
			v[1589] = -(v[1437] * v[1626] * v[2606]);
			v[1627] = v[1621] * v[2487] + v[2607] * v[482];
			v[1642] = v[1627] * v[2643];
			v[1629] = v[1623] * v[2487] + v[2607] * v[483];
			v[1646] = 2e0 * v[1629] * v[495];
			v[1630] = v[1624] * v[2487] + v[2607] * v[484];
			v[1643] = v[1630] * v[2640];
			v[1592] = v[1635] * v[2580] + v[2608] * v[487];
			v[1591] = v[1625] * v[1633] * v[491];
			v[1590] = v[1625] * v[2580] * v[492] + v[2608] * v[493];
			v[1616] = v[1635] * v[1636] * v[487];
			v[1617] = v[1436] * v[1625] * v[1631] * v[1636] * v[2773];
			v[1588] = v[1624] * v[2581] + v[1431] * v[2607];
			v[1587] = v[1623] * v[2581] + v[1433] * v[2607];
			v[1586] = v[1621] * v[2581] + v[1434] * v[2607];
			v[1638] = (v[1629] * v[494] + v[1627] * v[495]) / 2e0;
			v[1639] = v[1642] + v[1646];
			v[1640] = v[1639] + v[1643];
			v[1641] = (v[1630] * v[494] + v[1627] * v[496]) / 2e0;
			v[1644] = v[1642] + v[1643];
			v[1645] = (v[1630] * v[495] + v[1629] * v[496]) / 2e0;
			v[1647] = v[1643] + v[1646];
			v[1595] = (v[1417] * v[1627] + v[1404] * v[1629] + 4e0 * v[1630] * v[1648]) / 2e0;
			v[1594] = (v[1424] * v[1627] + v[1404] * v[1630] + v[1629] * v[2609]) / 2e0;
			v[1593] = (v[1424] * v[1629] + v[1417] * v[1630] + v[1627] * v[2610]) / 2e0;
			v[1651] = v[1640] * v[2639];
			v[1615] = 8e0 * v[1426] * v[1640] * v[2774];
			v[1614] = v[1651] * v[2611];
			v[1652] = v[1630] + v[1638];
			v[1653] = v[1630] - v[1638];
			v[1613] = v[1350] * v[1651];
			v[1654] = -v[1629] + v[1641];
			v[1655] = v[1629] + v[1641];
			v[1612] = v[1351] * v[1651];
			v[1600] = v[1411] * v[1651] + v[1652] * v[497];
			v[1611] = v[1352] * v[1651];
			v[1601] = (v[1407] * v[1651] - v[1644] * v[497]) / 2e0;
			v[1610] = v[1651] * v[2612];
			v[1656] = v[1627] + v[1645];
			v[1657] = -v[1627] + v[1645];
			v[1609] = v[1354] * v[1651];
			v[1603] = v[1398] * v[1651] + v[1654] * v[497];
			v[1608] = v[1355] * v[1651];
			v[1604] = v[1394] * v[1651] + v[1656] * v[497];
			v[1607] = v[1356] * v[1651];
			v[1605] = (v[1390] * v[1651] - v[1639] * v[497]) / 2e0;
			v[1606] = v[1651] * v[2613];
			v[1602] = v[1402] * v[1651] + v[1657] * v[497];
			v[1599] = v[1415] * v[1651] + v[1655] * v[497];
			v[1598] = v[1420] * v[1651] - v[1653] * v[497];
			v[1597] = (v[1425] * v[1651] - v[1647] * v[497]) / 2e0;
			v[1596] = v[1352] * v[1652] - v[1350] * v[1653] + v[1355] * v[1654] + v[1351] * v[1655] + v[1356] * v[1656]
				+ v[1354] * v[1657] - v[1647] * v[2611] - v[1644] * v[2612] - v[1639] * v[2613];
		}
		else {
		};
		v[1658] = 0e0;
		v[1659] = 0e0;
		v[1660] = 0e0;
		v[1661] = 0e0;
		v[1662] = 0e0;
		v[1663] = 0e0;
		b1664 = b23;
		if (b1664) {
			v[2616] = -(v[1585] * v[517]);
			v[2615] = -(v[1581] * v[516]);
			v[2614] = -(v[1577] * v[515]);
			v[1698] = v[1343] * v[1585];
			v[1671] = i2587 * v[160] + i2588 * v[167] + i2589 * v[174] - i2590 * v[262] - i2591 * v[269] - i2592 * v[276]
				+ v[1605] * v[7];
			v[1605] = 0e0;
			v[1672] = v[1671] + v[1604] * v[6];
			v[1604] = 0e0;
			v[1673] = v[1672] + v[1603] * v[5];
			v[2619] = -(v[1673] * v[360]);
			v[1603] = 0e0;
			v[1680] = i2593 * v[160] + i2594 * v[167] + i2595 * v[174] - i2596 * v[262] - i2597 * v[269] - i2598 * v[276]
				+ v[1602] * v[7];
			v[1602] = 0e0;
			v[1681] = v[1680] + v[1601] * v[6];
			v[1601] = 0e0;
			v[1682] = v[1681] + v[1600] * v[5];
			v[2621] = -(v[1682] * v[359]);
			v[1600] = 0e0;
			v[1689] = i2599 * v[160] + i2600 * v[167] + i2601 * v[174] - i2602 * v[262] - i2603 * v[269] - i2604 * v[276]
				+ v[1599] * v[7];
			v[1599] = 0e0;
			v[1690] = v[1689] + v[1598] * v[6];
			v[1598] = 0e0;
			v[1691] = v[1690] + v[1597] * v[5];
			v[2622] = -(v[1691] * v[358]);
			v[1597] = 0e0;
			v[1692] = v[1345] * v[1577];
			v[1658] = v[1577] * v[1702];
			v[1584] = v[1584] + v[1343] * v[2614];
			v[1580] = v[1580] + v[1344] * v[2614];
			v[1577] = 0e0;
			v[1694] = v[1344] * v[1581];
			v[1695] = v[1692] + v[1694];
			v[1659] = v[1581] * v[1703];
			v[1584] = v[1584] + v[1343] * v[2615];
			v[1576] = v[1576] + v[1345] * v[2615];
			v[1581] = 0e0;
			v[1697] = v[1694] + v[1698];
			v[1699] = v[1692] + v[1698];
			v[1660] = v[1585] * v[1704];
			v[1580] = v[1580] + v[1344] * v[2616];
			v[1576] = v[1576] + v[1345] * v[2616];
			v[1585] = 0e0;
			v[1663] = v[1343] * v[1673];
			v[1662] = v[1344] * v[1682];
			v[1661] = v[1345] * v[1691];
			v[1658] = v[1658] - (v[1697] + v[1345] * v[2617]) * v[358];
			v[1576] = v[1576] + v[1691] * v[1702] + v[1345] * (v[2619] + v[2621]) - v[1697] * v[515];
			v[1659] = v[1659] - (v[1699] + v[1344] * v[2618]) * v[359];
			v[1580] = v[1580] + v[1682] * v[1703] + v[1344] * (v[2619] + v[2622]) - v[1699] * v[516];
			v[1660] = v[1660] - (v[1695] + v[1343] * v[2620]) * v[360];
			v[1584] = v[1584] + v[1673] * v[1704] + v[1343] * (v[2621] + v[2622]) - v[1695] * v[517];
		}
		else {
		};
		v[1705] = v[14] * (v[1217] * (i2587 * v[544] + i2588 * v[547] + i2589 * v[552]) + v[1216] * (i2593 * v[544] + i2594 * v[547]
			+ i2595 * v[552]) + v[1215] * (i2599 * v[544] + i2600 * v[547] + i2601 * v[552]));
		v[1706] = v[14] * (v[1217] * (i2590 * v[562] + i2591 * v[565] + i2592 * v[570]) + v[1216] * (i2596 * v[562] + i2597 * v[565]
			+ i2598 * v[570]) + v[1215] * (i2602 * v[562] + i2603 * v[565] + i2604 * v[570]));
		v[1725] = -(i2599 * v[1707]) - i2593 * v[1708] - i2587 * v[1709] - i2600 * v[1710] - i2594 * v[1711] - i2588 * v[1712]
			- i2601 * v[1713] - i2595 * v[1714] - i2589 * v[1715] - i2602 * v[1716] - i2596 * v[1717] - i2590 * v[1718] - i2603 * v[1719]
			- i2597 * v[1720] - i2591 * v[1721] - i2604 * v[1722] - i2598 * v[1723] - i2592 * v[1724];
		v[1782] = v[1725];
		v[1744] = -(i2599 * v[1726]) - i2593 * v[1727] - i2587 * v[1728] - i2600 * v[1729] - i2594 * v[1730] - i2588 * v[1731]
			- i2601 * v[1732] - i2595 * v[1733] - i2589 * v[1734] - i2602 * v[1735] - i2596 * v[1736] - i2590 * v[1737] - i2603 * v[1738]
			- i2597 * v[1739] - i2591 * v[1740] - i2604 * v[1741] - i2598 * v[1742] - i2592 * v[1743];
		v[1780] = v[1744];
		v[1763] = -(i2599 * v[1745]) - i2593 * v[1746] - i2587 * v[1747] - i2600 * v[1748] - i2594 * v[1749] - i2588 * v[1750]
			- i2601 * v[1751] - i2595 * v[1752] - i2589 * v[1753] - i2602 * v[1754] - i2596 * v[1755] - i2590 * v[1756] - i2603 * v[1757]
			- i2597 * v[1758] - i2591 * v[1759] - i2604 * v[1760] - i2598 * v[1761] - i2592 * v[1762];
		v[1778] = v[1763];
		v[1764] = v[2623] * v[600];
		v[1765] = v[2623] * v[601];
		v[1766] = v[2623] * v[602];
		v[1767] = -(v[2624] * v[600]);
		v[1768] = -(v[2624] * v[601]);
		v[1769] = -(v[2624] * v[602]);
		v[1770] = 0e0;
		v[1771] = 0e0;
		v[1772] = 0e0;
		v[1773] = 0e0;
		v[1774] = 0e0;
		v[1775] = 0e0;
		b1776 = b22;
		if (b1776) {
			b1777 = b598;
			if (b1777) {
				v[1775] = v[1763];
				v[1763] = 0e0;
				v[1774] = v[1744];
				v[1744] = 0e0;
				v[1773] = v[1725];
				v[1725] = 0e0;
			}
			else {
				v[1790] = v[1782] * v[2495];
				v[1788] = v[1780] * v[2495];
				v[1786] = v[1778] * v[2495];
				v[1763] = 0e0;
				v[1744] = 0e0;
				v[2626] = (v[12] * (v[1782] * v[610] + v[1780] * v[611] + v[1778] * v[612])) / sqrt(v[1791]);
				v[1725] = 0e0;
				v[2625] = ((v[1790] * v[594] + v[1788] * v[595] + v[1786] * v[596]) * v[608]) / sqrt(v[1787]);
				v[1775] = v[2625] * v[596] + v[1786] * v[609];
				v[1774] = v[2625] * v[595] + v[1788] * v[609];
				v[1773] = v[2625] * v[594] + v[1790] * v[609];
				v[1772] = v[2626] * v[587];
				v[1771] = v[2626] * v[586];
				v[1770] = v[2626] * v[585];
			};
		}
		else {
			b1793 = b632;
			if (b1793) {
				v[1775] = v[1778];
				v[1763] = 0e0;
				v[1774] = v[1780];
				v[1744] = 0e0;
				v[1773] = v[1782];
				v[1725] = 0e0;
			}
			else {
				v[1799] = v[12] * v[1782] * v[643];
				v[1797] = v[12] * v[1780] * v[643];
				v[1796] = v[12] * v[1778] * v[643];
				v[2628] = (v[12] * (v[1782] * v[640] + v[1780] * v[641] + v[1778] * v[642])) / sqrt(v[1791]);
				v[2627] = ((v[1799] * v[594] + v[1797] * v[595] + v[1796] * v[596]) * v[638]) / sqrt(v[1787]);
				v[1775] = v[2627] * v[596] + v[1796] * v[639];
				v[1774] = v[2627] * v[595] + v[1797] * v[639];
				v[1773] = v[2627] * v[594] + v[1799] * v[639];
				v[1772] = v[2628] * v[587];
				v[1771] = v[2628] * v[586];
				v[1770] = v[2628] * v[585];
			};
		};
		v[2632] = v[13] * v[1770];
		v[2002] = -(v[2632] * v[520]);
		v[2631] = v[13] * v[1771];
		v[1968] = -(v[2631] * v[528]);
		v[2630] = v[13] * v[1772];
		v[1805] = -(v[14] * (v[1217] * v[1559] + v[1775] * v[340]));
		v[1806] = -(v[14] * (v[1217] * v[1558] + v[1775] * v[339]));
		v[1807] = -(v[14] * (v[1217] * v[1557] + v[1775] * v[338]));
		v[1808] = -(v[14] * (v[1217] * v[1556] + v[1775] * v[337]));
		v[1809] = -(v[14] * (v[1217] * v[1555] + v[1775] * v[336]));
		v[1810] = -(v[14] * (v[1217] * v[1554] + v[1775] * v[335]));
		v[1811] = -(v[14] * (v[1217] * v[1553] + v[1775] * v[334]));
		v[1812] = -(v[14] * (v[1217] * v[1552] + v[1775] * v[333]));
		v[1813] = -(v[14] * (v[1217] * v[1551] + v[1775] * v[332]));
		v[1814] = -(v[14] * (v[1217] * v[1541] + v[1775] * v[331]));
		v[1815] = -(v[14] * (v[1217] * v[1540] + v[1775] * v[330]));
		v[1816] = -(v[14] * (v[1217] * v[1539] + v[1775] * v[329]));
		v[1817] = -(v[14] * (v[1217] * v[1538] + v[1775] * v[328]));
		v[1818] = -(v[14] * (v[1217] * v[1537] + v[1775] * v[327]));
		v[1819] = -(v[14] * (v[1217] * v[1536] + v[1775] * v[326]));
		v[1820] = -(v[14] * (v[1217] * v[1535] + v[1775] * v[325]));
		v[1821] = -(v[14] * (v[1217] * v[1534] + v[1775] * v[324]));
		v[1822] = -(v[14] * (v[1217] * v[1533] + v[1775] * v[323]));
		v[1842] = -(v[14] * (v[1216] * v[1559] + v[1774] * v[340]));
		v[1843] = -(v[14] * (v[1216] * v[1558] + v[1774] * v[339]));
		v[1844] = -(v[14] * (v[1216] * v[1557] + v[1774] * v[338]));
		v[1845] = -(v[14] * (v[1216] * v[1556] + v[1774] * v[337]));
		v[1846] = -(v[14] * (v[1216] * v[1555] + v[1774] * v[336]));
		v[1847] = -(v[14] * (v[1216] * v[1554] + v[1774] * v[335]));
		v[1848] = -(v[14] * (v[1216] * v[1553] + v[1774] * v[334]));
		v[1849] = -(v[14] * (v[1216] * v[1552] + v[1774] * v[333]));
		v[1850] = -(v[14] * (v[1216] * v[1551] + v[1774] * v[332]));
		v[1851] = -(v[14] * (v[1216] * v[1541] + v[1774] * v[331]));
		v[1852] = -(v[14] * (v[1216] * v[1540] + v[1774] * v[330]));
		v[1853] = -(v[14] * (v[1216] * v[1539] + v[1774] * v[329]));
		v[1854] = -(v[14] * (v[1216] * v[1538] + v[1774] * v[328]));
		v[1855] = -(v[14] * (v[1216] * v[1537] + v[1774] * v[327]));
		v[1856] = -(v[14] * (v[1216] * v[1536] + v[1774] * v[326]));
		v[1857] = -(v[14] * (v[1216] * v[1535] + v[1774] * v[325]));
		v[1858] = -(v[14] * (v[1216] * v[1534] + v[1774] * v[324]));
		v[1859] = -(v[14] * (v[1216] * v[1533] + v[1774] * v[323]));
		v[1879] = -(v[14] * (v[1215] * v[1559] + v[1773] * v[340]));
		v[1880] = -(v[14] * (v[1215] * v[1558] + v[1773] * v[339]));
		v[1881] = -(v[14] * (v[1215] * v[1557] + v[1773] * v[338]));
		v[1882] = -(v[14] * (v[1215] * v[1556] + v[1773] * v[337]));
		v[1883] = -(v[14] * (v[1215] * v[1555] + v[1773] * v[336]));
		v[1884] = -(v[14] * (v[1215] * v[1554] + v[1773] * v[335]));
		v[1885] = -(v[14] * (v[1215] * v[1553] + v[1773] * v[334]));
		v[1886] = -(v[14] * (v[1215] * v[1552] + v[1773] * v[333]));
		v[1887] = -(v[14] * (v[1215] * v[1551] + v[1773] * v[332]));
		v[1888] = -(v[14] * (v[1215] * v[1541] + v[1773] * v[331]));
		v[1889] = -(v[14] * (v[1215] * v[1540] + v[1773] * v[330]));
		v[1890] = -(v[14] * (v[1215] * v[1539] + v[1773] * v[329]));
		v[1891] = -(v[14] * (v[1215] * v[1538] + v[1773] * v[328]));
		v[1892] = -(v[14] * (v[1215] * v[1537] + v[1773] * v[327]));
		v[1893] = -(v[14] * (v[1215] * v[1536] + v[1773] * v[326]));
		v[1894] = -(v[14] * (v[1215] * v[1535] + v[1773] * v[325]));
		v[1895] = -(v[14] * (v[1215] * v[1534] + v[1773] * v[324]));
		v[1896] = -(v[14] * (v[1215] * v[1533] + v[1773] * v[323]));
		v[1915] = v[10] * v[1775];
		v[1916] = v[10] * v[1774];
		v[1917] = v[10] * v[1773];
		v[1918] = v[1772] * v[360];
		v[2629] = v[13] * v[1918];
		v[1944] = v[200] * v[2629];
		v[1940] = v[205] * v[2629];
		v[1936] = v[211] * v[2629];
		v[1932] = -(v[2629] * v[302]);
		v[1928] = -(v[2629] * v[307]);
		v[1923] = -(v[2629] * v[313]);
		v[1919] = v[1008] * v[2630];
		v[1584] = v[1584] + v[1772] * v[1920];
		v[1580] = v[1580] + v[1006] * v[2630];
		v[1576] = v[1576] + v[1005] * v[2630];
		v[1965] = v[1023] * v[2631];
		v[1580] = v[1580] + v[1031] * v[1771];
		v[1985] = v[1918] + v[1770] * v[358] + v[1771] * v[359];
		v[1986] = v[13] * (-v[1918] + v[1985]);
		v[2650] = -(v[1986] * v[360]);
		v[2651] = v[2650] - v[2630] * v[536];
		v[1992] = v[13] * (v[1771] * v[323] + v[1770] * v[324]);
		v[1993] = v[13] * (v[1771] * v[326] + v[1770] * v[327]);
		v[1994] = v[13] * (v[1771] * v[329] + v[1770] * v[330]);
		v[1995] = v[13] * (v[1771] * v[332] + v[1770] * v[333]);
		v[1996] = v[13] * (v[1771] * v[335] + v[1770] * v[336]);
		v[1997] = v[13] * (v[1771] * v[338] + v[1770] * v[339]);
		v[2647] = v[1992] * v[200] + v[1993] * v[205] + v[1994] * v[211] - v[1995] * v[302] - v[1996] * v[307] - v[1997] * v[313];
		v[1998] = v[1022] * v[2632];
		v[1576] = v[1576] + v[1031] * v[1770];
		v[2018] = 0e0;
		v[2019] = 0e0;
		v[2020] = 0e0;
		v[2021] = 0e0;
		v[2022] = 0e0;
		v[2023] = 0e0;
		v[2024] = 0e0;
		v[2025] = 0e0;
		v[2026] = 0e0;
		b2027 = b23;
		if (b2027) {
			v[2635] = v[1917] * v[358];
			v[2634] = v[1916] * v[359];
			v[2636] = v[2634] + v[2635];
			v[2633] = v[1915] * v[360];
			v[2638] = v[2633] + v[2634];
			v[2637] = v[2633] + v[2635];
			v[1663] = v[1663] + v[1915] * v[517];
			v[1662] = v[1662] + v[1916] * v[516];
			v[1658] = v[1658] + v[1361] * v[1917] - v[2638] * v[358];
			v[1659] = v[1659] + v[1363] * v[1916] - v[2637] * v[359];
			v[1660] = v[1660] + v[1365] * v[1915] - v[2636] * v[360];
			v[1661] = v[1661] + v[1917] * v[515];
			v[1584] = v[1584] + v[1915] * v[2576] - v[2636] * v[517];
			v[1580] = v[1580] + v[1916] * v[2577] - v[2637] * v[516];
			v[1576] = v[1576] - v[1917] * v[2578] - v[2638] * v[515];
			v[2018] = v[1658] * v[5];
			v[2019] = v[1658] * v[6];
			v[2020] = v[1658] * v[7];
			v[2028] = -(v[1658] * v[276]);
			v[2029] = -(v[1658] * v[269]);
			v[2030] = -(v[1658] * v[262]);
			v[2031] = v[1658] * v[174];
			v[2032] = v[1658] * v[167];
			v[2033] = v[160] * v[1658];
			v[2021] = v[1659] * v[5];
			v[2022] = v[1659] * v[6];
			v[2023] = v[1659] * v[7];
			v[2034] = -(v[1659] * v[276]);
			v[2035] = -(v[1659] * v[269]);
			v[2036] = -(v[1659] * v[262]);
			v[2037] = v[1659] * v[174];
			v[2038] = v[1659] * v[167];
			v[2039] = v[160] * v[1659];
			v[2024] = v[1660] * v[5];
			v[2025] = v[1660] * v[6];
			v[2026] = v[1660] * v[7];
			v[2040] = -(v[1660] * v[276]);
			v[2041] = -(v[1660] * v[269]);
			v[2042] = -(v[1660] * v[262]);
			v[2043] = v[1660] * v[174];
			v[2044] = v[1660] * v[167];
			v[2045] = v[160] * v[1660];
			v[1998] = -v[1661] + v[1998];
			v[1965] = -v[1662] + v[1965];
			v[1919] = -v[1663] + v[1919];
		}
		else {
			v[2033] = 0e0;
			v[2039] = 0e0;
			v[2045] = 0e0;
			v[2032] = 0e0;
			v[2038] = 0e0;
			v[2044] = 0e0;
			v[2031] = 0e0;
			v[2037] = 0e0;
			v[2043] = 0e0;
			v[2030] = 0e0;
			v[2036] = 0e0;
			v[2042] = 0e0;
			v[2029] = 0e0;
			v[2035] = 0e0;
			v[2041] = 0e0;
			v[2028] = 0e0;
			v[2034] = 0e0;
			v[2040] = 0e0;
		};
		b2046 = b480;
		if (b2046) {
			v[1596] = v[1596] + (v[1390] * v[2026]) / 2e0;
			v[1606] = v[1606] + v[2026] * v[2489];
			v[1596] = v[1596] + v[1394] * v[2025];
			v[1607] = v[1607] + v[2025] * v[497];
			v[1596] = v[1596] + v[1398] * v[2024];
			v[1608] = v[1608] + v[2024] * v[497];
			v[1596] = v[1596] + v[1402] * v[2023];
			v[1609] = v[1609] + v[2023] * v[497];
			v[1596] = v[1596] + (v[1407] * v[2022]) / 2e0;
			v[1610] = v[1610] + v[2022] * v[2489];
			v[1596] = v[1596] + v[1411] * v[2021];
			v[1611] = v[1611] + v[2021] * v[497];
			v[1596] = v[1596] + v[1415] * v[2020];
			v[1612] = v[1612] + v[2020] * v[497];
			v[1596] = v[1596] + v[1420] * v[2019];
			v[1613] = v[1613] + v[2019] * v[497];
			v[1596] = v[1596] + (v[1425] * v[2018]) / 2e0;
			v[1614] = v[1614] + v[2018] * v[2489];
			v[1615] = v[1615] + v[1596] * v[2639];
			v[2644] = -v[1610] + v[1615];
			v[1594] = v[1594] - v[1608];
			v[2056] = v[1608] + v[1612];
			v[1594] = v[1594] + v[1612];
			v[1593] = v[1593] + v[1607] + (v[2056] * v[496]) / 2e0;
			v[2058] = v[1607] + v[1609];
			v[1593] = v[1593] - v[1609];
			v[1595] = v[1595] + v[1611] + v[2058] * v[2488] + v[2056] * v[2579] + v[2640] * (-v[1614] + v[2644]);
			v[1595] = v[1595] - v[1613];
			v[2641] = v[1595] * v[484];
			v[2060] = v[1611] + v[1613];
			v[1592] = v[1592] + v[2641] * v[487];
			v[1590] = v[1590] + v[2641] * v[493];
			v[1588] = v[1588] + v[1595] * v[2487];
			v[1594] = v[1594] + (v[2060] * v[494] - 4e0 * (v[1606] + v[1614] - v[1615]) * v[495] + v[2058] * v[496]) / 2e0;
			v[2642] = v[1594] * v[483];
			v[1592] = v[1592] + v[2642] * v[487];
			v[1590] = v[1590] + v[2642] * v[493];
			v[1587] = v[1587] + v[1594] * v[2487];
			v[1593] = v[1593] + v[2060] * v[2488] + v[2643] * (-v[1606] + v[2644]);
			v[2645] = v[1593] * v[482];
			v[1592] = v[1592] + v[2645] * v[487];
			v[1590] = v[1590] + v[2645] * v[493];
			v[1586] = v[1586] + v[1593] * v[2487];
			v[1591] = v[1591] + v[1592] * v[492];
			v[1589] = v[1589] + v[1591];
			v[1616] = v[1616] + 2e0 * v[1590] * v[1631];
			v[1617] = v[1617] + (v[1616] * v[1632]) / 2e0;
			v[1589] = v[1589] + v[1617];
			v[2646] = v[1589] / v[485];
			v[1588] = v[1588] + v[2646] * v[484];
			v[1587] = v[1587] + v[2646] * v[483];
			v[1586] = v[1586] + v[2646] * v[482];
			v[1580] = v[1580] + v[1588] * v[367];
			v[1576] = v[1576] - v[1588] * v[368];
			v[1584] = v[1584] - v[1587] * v[367];
			v[1576] = v[1576] + v[1587] * v[369];
			v[1584] = v[1584] + v[1586] * v[368];
			v[1580] = v[1580] - v[1586] * v[369];
		}
		else {
		};
		v[2066] = -(v[1822] * v[462]) - v[1821] * v[463] - v[1820] * v[464] - v[1819] * v[465] - v[1818] * v[466]
			- v[1817] * v[467] - v[1816] * v[468] - v[1815] * v[469] - v[1814] * v[470] - v[1813] * v[471] - v[1812] * v[472]
			- v[1811] * v[473] - v[1810] * v[474] - v[1809] * v[475] - v[1808] * v[476] - v[1807] * v[477] - v[1806] * v[478]
			- v[1805] * v[479];
		v[2067] = v[1822] * v[444] + v[1821] * v[445] + v[1820] * v[446] + v[1819] * v[447] + v[1818] * v[448] + v[1817] * v[449]
			+ v[1816] * v[450] + v[1815] * v[451] + v[1814] * v[452] + v[1813] * v[453] + v[1812] * v[454] + v[1811] * v[455]
			+ v[1810] * v[456] + v[1809] * v[457] + v[1808] * v[458] + v[1807] * v[459] + v[1806] * v[460] + v[1805] * v[461];
		v[2068] = -(v[1859] * v[462]) - v[1858] * v[463] - v[1857] * v[464] - v[1856] * v[465] - v[1855] * v[466]
			- v[1854] * v[467] - v[1853] * v[468] - v[1852] * v[469] - v[1851] * v[470] - v[1850] * v[471] - v[1849] * v[472]
			- v[1848] * v[473] - v[1847] * v[474] - v[1846] * v[475] - v[1845] * v[476] - v[1844] * v[477] - v[1843] * v[478]
			- v[1842] * v[479];
		v[2069] = v[1859] * v[444] + v[1858] * v[445] + v[1857] * v[446] + v[1856] * v[447] + v[1855] * v[448] + v[1854] * v[449]
			+ v[1853] * v[450] + v[1852] * v[451] + v[1851] * v[452] + v[1850] * v[453] + v[1849] * v[454] + v[1848] * v[455]
			+ v[1847] * v[456] + v[1846] * v[457] + v[1845] * v[458] + v[1844] * v[459] + v[1843] * v[460] + v[1842] * v[461];
		v[2070] = -(v[1896] * v[462]) - v[1895] * v[463] - v[1894] * v[464] - v[1893] * v[465] - v[1892] * v[466]
			- v[1891] * v[467] - v[1890] * v[468] - v[1889] * v[469] - v[1888] * v[470] - v[1887] * v[471] - v[1886] * v[472]
			- v[1885] * v[473] - v[1884] * v[474] - v[1883] * v[475] - v[1882] * v[476] - v[1881] * v[477] - v[1880] * v[478]
			- v[1879] * v[479];
		v[2071] = v[1896] * v[444] + v[1895] * v[445] + v[1894] * v[446] + v[1893] * v[447] + v[1892] * v[448] + v[1891] * v[449]
			+ v[1890] * v[450] + v[1889] * v[451] + v[1888] * v[452] + v[1887] * v[453] + v[1886] * v[454] + v[1885] * v[455]
			+ v[1884] * v[456] + v[1883] * v[457] + v[1882] * v[458] + v[1881] * v[459] + v[1880] * v[460] + v[1879] * v[461];
		v[1584] = v[1584] + v[1008] * v[1986] + v[1919] * v[2582];
		v[1580] = v[1580] + v[1965] * v[2583] + v[2647] * v[358];
		v[1576] = v[1576] + v[1998] * v[2584] + v[2647] * v[359];
		v[2081] = -(v[1034] * v[1448] * v[2648]) + v[1985] * v[2649] + v[1446] * v[1572] * v[355] + (v[1348] * v[1562]
			+ v[1347] * v[1566] + v[1346] * v[1570] + v[1576] * v[341] + v[1580] * v[342] + v[1584] * v[343]) * v[356];
		v[2083] = v[1348] * v[1573] + (v[1448] * v[1562] + v[2081] * v[343]) / v[351] + v[1584] * v[357];
		v[2085] = v[1347] * v[1573] + (v[1448] * v[1566] + v[2081] * v[342]) / v[351] + v[1580] * v[357];
		v[2087] = v[1346] * v[1573] + (v[1448] * v[1570] + v[2081] * v[341]) / v[351] + v[1576] * v[357];
		v[2088] = -(i2592 * v[1449]) - i2598 * v[1450] - i2604 * v[1451] - v[2087] * v[227] - v[2085] * v[228] - v[2083] * v[229]
			+ v[1997] * v[2479] + v[2002] * v[338] + v[1968] * v[339] + v[2650] * v[340] + v[2630] * v[979];
		v[2089] = -(i2591 * v[1449]) - i2597 * v[1450] - i2603 * v[1451] - v[2087] * v[224] - v[2085] * v[225] - v[2083] * v[226]
			+ v[1996] * v[2479] + v[2002] * v[335] + v[1968] * v[336] + v[2650] * v[337] + v[2630] * v[960];
		v[2090] = -(i2590 * v[1449]) - i2596 * v[1450] - i2602 * v[1451] - v[2087] * v[221] - v[2085] * v[222] - v[2083] * v[223]
			+ v[1995] * v[2479] + v[2002] * v[332] + v[1968] * v[333] + v[2650] * v[334] + v[2630] * v[941];
		v[2091] = i2589 * v[1449] + i2595 * v[1450] + i2601 * v[1451] + v[127] * v[2083] + v[126] * v[2085] + v[125] * v[2087]
			- v[1994] * v[2479] - v[2002] * v[329] - v[1968] * v[330] - v[2650] * v[331] + v[2630] * v[922];
		v[2092] = i2588 * v[1449] + i2594 * v[1450] + i2600 * v[1451] + v[124] * v[2083] + v[123] * v[2085] + v[122] * v[2087]
			- v[1993] * v[2479] - v[2002] * v[326] - v[1968] * v[327] - v[2650] * v[328] + v[2630] * v[900];
		v[2093] = i2587 * v[1449] + i2593 * v[1450] + i2599 * v[1451] + v[121] * v[2083] + v[120] * v[2085] + v[119] * v[2087]
			- v[1992] * v[2479] - v[2002] * v[323] - v[1968] * v[324] - v[2650] * v[325] + v[2630] * v[869];
		v[2094] = i2591 * v[1261] + i2597 * v[1301] + i2603 * v[1341] + v[2070] * v[224] + v[2068] * v[225] + v[2066] * v[226];
		v[2095] = i2592 * v[1261] + i2598 * v[1301] + i2604 * v[1341] + v[2070] * v[227] + v[2068] * v[228] + v[2066] * v[229];
		v[2028] = v[1190] * v[1764] + v[2028] - v[2087] * v[313] + v[2070] * v[570];
		v[2029] = v[1186] * v[1764] + v[2029] - v[2087] * v[307] + v[2070] * v[565];
		v[2030] = v[1183] * v[1764] + v[2030] - v[2087] * v[302] + v[2070] * v[562];
		v[2034] = v[1190] * v[1765] + v[2034] - v[2085] * v[313] + v[2068] * v[570];
		v[2035] = v[1186] * v[1765] + v[2035] - v[2085] * v[307] + v[2068] * v[565];
		v[2036] = v[1183] * v[1765] + v[2036] - v[2085] * v[302] + v[2068] * v[562];
		v[2096] = v[1764] * v[224] + v[1765] * v[225] + v[1766] * v[226];
		v[2040] = v[1190] * v[1766] + v[2040] - v[2083] * v[313] + v[2066] * v[570];
		v[2041] = v[1186] * v[1766] + v[2041] - v[2083] * v[307] + v[2066] * v[565];
		v[2042] = v[1183] * v[1766] + v[2042] - v[2083] * v[302] + v[2066] * v[562];
		v[2097] = v[1764] * v[227] + v[1765] * v[228] + v[1766] * v[229];
		v[2098] = 0e0;
		b2099 = b273;
		if (b2099) {
			v[2098] = v[2088] * v[568];
		}
		else {
		};
		v[2100] = 0e0;
		v[2101] = 0e0;
		b2102 = b266;
		if (b2102) {
			v[2100] = v[2088] * v[2465];
			v[2101] = -v[2095] - v[2097];
			v[2098] = v[2098] - v[2088] * v[400] + v[271] * (v[1179] * v[2097] + v[2095] * v[557]);
		}
		else {
		};
		b2103 = b266;
		if (b2103) {
			v[2101] = v[2094] + v[2096] + v[2101];
			v[2098] = v[2098] + v[2089] * v[400] - v[271] * (v[1179] * v[2096] + v[2094] * v[557]);
			v[2100] = v[2100] + v[271] * (-v[2101] + v[2089] * v[397]);
		}
		else {
		};
		v[2104] = 0e0;
		v[2105] = 0e0;
		b2106 = b259;
		if (b2106) {
			v[2104] = v[2089] * v[2463];
			v[2105] = -v[2094] - v[2096];
			v[2098] = v[2098] - v[2089] * v[394] + v[264] * (v[1178] * v[2096] + v[2094] * v[556]);
		}
		else {
		};
		v[2107] = i2590 * v[1261] + i2596 * v[1301] + i2602 * v[1341] + v[2070] * v[221] + v[2068] * v[222] + v[2066] * v[223];
		v[2108] = v[1764] * v[221] + v[1765] * v[222] + v[1766] * v[223];
		b2109 = b259;
		if (b2109) {
			v[2105] = v[2105] + v[2107] + v[2108];
			v[2104] = v[2104] + v[264] * (-v[2105] + v[2090] * v[730]);
			v[2098] = v[2098] + v[2090] * v[394] - v[264] * (v[1178] * v[2108] + v[2107] * v[556]);
		}
		else {
		};
		b2110 = b256;
		if (b2110) {
			v[2098] = v[2098] + v[2090] * v[558];
		}
		else {
		};
		v[2111] = v[2100];
		b2112 = b239;
		if (b2112) {
			v[2098] = v[2098] + v[2111] * v[244];
		}
		else {
		};
		v[2113] = v[2104];
		b2114 = b239;
		if (b2114) {
			v[2098] = v[2098] - v[2113] * v[244];
		}
		else {
		};
		v[2124] = i2588 * v[1262] + i2594 * v[1303] + i2600 * v[1342] + v[124] * v[2067] + v[123] * v[2069] + v[122] * v[2071];
		v[2125] = i2589 * v[1262] + i2595 * v[1303] + i2601 * v[1342] + v[127] * v[2067] + v[126] * v[2069] + v[125] * v[2071];
		v[2031] = v[1174] * v[1767] + v[2031] + v[2087] * v[211] + v[2071] * v[552];
		v[2032] = v[1170] * v[1767] + v[2032] + v[205] * v[2087] + v[2071] * v[547];
		v[2033] = v[1167] * v[1767] + v[2033] + v[200] * v[2087] + v[2071] * v[544];
		v[2037] = v[1174] * v[1768] + v[2037] + v[2085] * v[211] + v[2069] * v[552];
		v[2038] = v[1170] * v[1768] + v[2038] + v[205] * v[2085] + v[2069] * v[547];
		v[2039] = v[1167] * v[1768] + v[2039] + v[200] * v[2085] + v[2069] * v[544];
		v[2126] = v[122] * v[1767] + v[123] * v[1768] + v[124] * v[1769];
		v[2043] = v[1174] * v[1769] + v[2043] + v[2083] * v[211] + v[2067] * v[552];
		v[2044] = v[1170] * v[1769] + v[2044] + v[205] * v[2083] + v[2067] * v[547];
		v[2045] = v[1167] * v[1769] + v[2045] + v[200] * v[2083] + v[2067] * v[544];
		v[2127] = v[125] * v[1767] + v[126] * v[1768] + v[127] * v[1769];
		v[2128] = 0e0;
		b2129 = b171;
		if (b2129) {
			v[2128] = v[2091] * v[550];
		}
		else {
		};
		v[2130] = 0e0;
		v[2131] = 0e0;
		b2132 = b164;
		if (b2132) {
			v[2130] = v[2091] * v[2459];
			v[2131] = -v[2125] - v[2127];
			v[2128] = v[2128] - v[2091] * v[381] + v[169] * (v[1163] * v[2127] + v[2125] * v[539]);
		}
		else {
		};
		b2133 = b164;
		if (b2133) {
			v[2131] = v[2124] + v[2126] + v[2131];
			v[2128] = v[2128] + v[2092] * v[381] - v[169] * (v[1163] * v[2126] + v[2124] * v[539]);
			v[2130] = v[2130] + v[169] * (-v[2131] + v[2092] * v[378]);
		}
		else {
		};
		v[2134] = 0e0;
		v[2135] = 0e0;
		b2136 = b157;
		if (b2136) {
			v[2134] = v[2092] * v[2457];
			v[2135] = -v[2124] - v[2126];
			v[2128] = v[2128] - v[2092] * v[375] + v[162] * (v[1162] * v[2126] + v[2124] * v[538]);
		}
		else {
		};
		v[2137] = i2587 * v[1262] + i2593 * v[1303] + i2599 * v[1342] + v[121] * v[2067] + v[120] * v[2069] + v[119] * v[2071];
		v[2138] = v[119] * v[1767] + v[120] * v[1768] + v[121] * v[1769];
		b2139 = b157;
		if (b2139) {
			v[2135] = v[2135] + v[2137] + v[2138];
			v[2134] = v[2134] + v[162] * (-v[2135] + v[2093] * v[745]);
			v[2128] = v[2128] + v[2093] * v[375] - v[162] * (v[1162] * v[2138] + v[2137] * v[538]);
		}
		else {
		};
		b2140 = b154;
		if (b2140) {
			v[2128] = v[2128] + v[2093] * v[540];
		}
		else {
		};
		v[2141] = v[2130];
		b2142 = b137;
		if (b2142) {
			v[2128] = v[2128] + v[142] * v[2141];
		}
		else {
		};
		v[2143] = v[2134];
		b2144 = b137;
		if (b2144) {
			v[2128] = v[2128] - v[142] * v[2143];
		}
		else {
		};
		v[4446] = v[2033] + v[2128] * v[444] + v[2098] * v[462] + v[15] * (v[14] * (-(v[1339] * v[1773]) - v[1298] * v[1774]
			- v[1259] * v[1775]) - v[200] * v[2002] + v[1944] * v[358] - v[1705] * v[444] + v[1706] * v[462] + v[2631] * v[522]);
		v[4447] = v[2039] + v[2128] * v[445] + v[2098] * v[463] + v[15] * (v[14] * (-(v[1337] * v[1773]) - v[1296] * v[1774]
			- v[1257] * v[1775]) - v[1968] * v[200] + v[1944] * v[359] - v[1705] * v[445] + v[1706] * v[463] + v[2632] * v[522]);
		v[4448] = v[2045] + v[2128] * v[446] + v[2098] * v[464] + v[15] * (v[14] * (-(v[1335] * v[1773]) - v[1294] * v[1774]
			- v[1255] * v[1775]) - v[200] * v[2651] - v[1705] * v[446] + v[1706] * v[464]);
		v[4449] = v[2032] + v[2128] * v[447] + v[2098] * v[465] + v[15] * (v[14] * (-(v[1333] * v[1773]) - v[1292] * v[1774]
			- v[1253] * v[1775]) - v[2002] * v[205] + v[1940] * v[358] - v[1705] * v[447] + v[1706] * v[465] + v[2631] * v[523]);
		v[4450] = v[2038] + v[2128] * v[448] + v[2098] * v[466] + v[15] * (v[14] * (-(v[1331] * v[1773]) - v[1290] * v[1774]
			- v[1251] * v[1775]) - v[1968] * v[205] + v[1940] * v[359] - v[1705] * v[448] + v[1706] * v[466] + v[2632] * v[523]);
		v[4451] = v[2044] + v[2128] * v[449] + v[2098] * v[467] + v[15] * (v[14] * (-(v[1329] * v[1773]) - v[1288] * v[1774]
			- v[1249] * v[1775]) - v[205] * v[2651] - v[1705] * v[449] + v[1706] * v[467]);
		v[4452] = v[2031] + v[2128] * v[450] + v[2098] * v[468] + v[15] * (v[14] * (-(v[1327] * v[1773]) - v[1286] * v[1774]
			- v[1247] * v[1775]) - v[2002] * v[211] + v[1936] * v[358] - v[1705] * v[450] + v[1706] * v[468] + v[2631] * v[524]);
		v[4453] = v[2037] + v[2128] * v[451] + v[2098] * v[469] + v[15] * (v[14] * (-(v[1325] * v[1773]) - v[1284] * v[1774]
			- v[1245] * v[1775]) - v[1968] * v[211] + v[1936] * v[359] - v[1705] * v[451] + v[1706] * v[469] + v[2632] * v[524]);
		v[4454] = v[2043] + v[2128] * v[452] + v[2098] * v[470] + v[15] * (v[14] * (-(v[1323] * v[1773]) - v[1282] * v[1774]
			- v[1243] * v[1775]) - v[211] * v[2651] - v[1705] * v[452] + v[1706] * v[470]);
		v[4455] = v[2030] + v[2128] * v[453] + v[2098] * v[471] + v[15] * (v[14] * (-(v[1321] * v[1773]) - v[1280] * v[1774]
			- v[1241] * v[1775]) + v[2002] * v[302] + v[1932] * v[358] - v[1705] * v[453] + v[1706] * v[471] + v[2631] * v[525]);
		v[4456] = v[2036] + v[2128] * v[454] + v[2098] * v[472] + v[15] * (v[14] * (-(v[1319] * v[1773]) - v[1278] * v[1774]
			- v[1239] * v[1775]) + v[1968] * v[302] + v[1932] * v[359] - v[1705] * v[454] + v[1706] * v[472] + v[2632] * v[525]);
		v[4457] = v[2042] + v[2128] * v[455] + v[2098] * v[473] + v[15] * (v[14] * (-(v[1317] * v[1773]) - v[1276] * v[1774]
			- v[1237] * v[1775]) + v[2651] * v[302] - v[1705] * v[455] + v[1706] * v[473]);
		v[4458] = v[2029] + v[2128] * v[456] + v[2098] * v[474] + v[15] * (v[14] * (-(v[1315] * v[1773]) - v[1274] * v[1774]
			- v[1235] * v[1775]) + v[2002] * v[307] + v[1928] * v[358] - v[1705] * v[456] + v[1706] * v[474] + v[2631] * v[526]);
		v[4459] = v[2035] + v[2128] * v[457] + v[2098] * v[475] + v[15] * (v[14] * (-(v[1313] * v[1773]) - v[1272] * v[1774]
			- v[1233] * v[1775]) + v[1968] * v[307] + v[1928] * v[359] - v[1705] * v[457] + v[1706] * v[475] + v[2632] * v[526]);
		v[4460] = v[2041] + v[2128] * v[458] + v[2098] * v[476] + v[15] * (v[14] * (-(v[1311] * v[1773]) - v[1270] * v[1774]
			- v[1231] * v[1775]) + v[2651] * v[307] - v[1705] * v[458] + v[1706] * v[476]);
		v[4461] = v[2028] + v[2128] * v[459] + v[2098] * v[477] + v[15] * (v[14] * (-(v[1309] * v[1773]) - v[1268] * v[1774]
			- v[1229] * v[1775]) + v[2002] * v[313] + v[1923] * v[358] - v[1705] * v[459] + v[1706] * v[477] + v[2631] * v[527]);
		v[4462] = v[2034] + v[2128] * v[460] + v[2098] * v[478] + v[15] * (v[14] * (-(v[1307] * v[1773]) - v[1266] * v[1774]
			- v[1227] * v[1775]) + v[1968] * v[313] + v[1923] * v[359] - v[1705] * v[460] + v[1706] * v[478] + v[2632] * v[527]);
		v[4463] = v[2040] + v[2128] * v[461] + v[2098] * v[479] + v[15] * (v[14] * (-(v[1305] * v[1773]) - v[1264] * v[1774]
			- v[1225] * v[1775]) + v[2651] * v[313] - v[1705] * v[461] + v[1706] * v[479]);
		Rc[i1213 - 1] += v[4423 + i1213];
		for (i1513 = 1; i1513 <= 18; i1513++) {
			Kc[i1213 - 1][i1513 - 1] += v[4445 + i1513];
		};/* end for */
	};/* end for */
	v[2196] = -(v[227] * v[582]) - v[228] * v[583] - v[229] * v[584];
	v[2197] = -(v[224] * v[582]) - v[225] * v[583] - v[226] * v[584];
	v[2198] = -(v[221] * v[582]) - v[222] * v[583] - v[223] * v[584];
	v[2199] = v[125] * v[582] + v[126] * v[583] + v[127] * v[584];
	v[2200] = v[122] * v[582] + v[123] * v[583] + v[124] * v[584];
	v[2201] = v[119] * v[582] + v[120] * v[583] + v[121] * v[584];
	v[2202] = 0e0;
	b2203 = b273;
	if (b2203) {
		v[2202] = v[2196] * v[568];
	}
	else {
	};
	v[2204] = 0e0;
	b2205 = b266;
	if (b2205) {
		v[2204] = v[2196] * v[2465];
		v[2202] = v[2202] - v[2196] * v[400];
	}
	else {
	};
	b2206 = b266;
	if (b2206) {
		v[2204] = v[2204] + v[2197] * v[2464];
		v[2202] = v[2202] + v[2197] * v[400];
	}
	else {
	};
	v[2207] = 0e0;
	b2208 = b259;
	if (b2208) {
		v[2207] = v[2197] * v[2463];
		v[2202] = v[2202] - v[2197] * v[394];
	}
	else {
	};
	b2209 = b259;
	if (b2209) {
		v[2207] = v[2207] + v[2198] * v[2462];
		v[2202] = v[2202] + v[2198] * v[394];
	}
	else {
	};
	b2210 = b256;
	if (b2210) {
		v[2202] = v[2202] + v[2198] * v[558];
	}
	else {
	};
	v[2211] = v[2204];
	b2212 = b239;
	if (b2212) {
		v[2202] = v[2202] + v[2211] * v[244];
	}
	else {
	};
	v[2213] = v[2207];
	b2214 = b239;
	if (b2214) {
		v[2202] = v[2202] - v[2213] * v[244];
	}
	else {
	};
	v[2215] = 0e0;
	b2216 = b171;
	if (b2216) {
		v[2215] = v[2199] * v[550];
	}
	else {
	};
	v[2217] = 0e0;
	b2218 = b164;
	if (b2218) {
		v[2217] = v[2199] * v[2459];
		v[2215] = v[2215] - v[2199] * v[381];
	}
	else {
	};
	b2219 = b164;
	if (b2219) {
		v[2217] = v[2217] + v[2200] * v[2458];
		v[2215] = v[2215] + v[2200] * v[381];
	}
	else {
	};
	v[2220] = 0e0;
	b2221 = b157;
	if (b2221) {
		v[2220] = v[2200] * v[2457];
		v[2215] = v[2215] - v[2200] * v[375];
	}
	else {
	};
	b2222 = b157;
	if (b2222) {
		v[2220] = v[2220] + v[2201] * v[2456];
		v[2215] = v[2215] + v[2201] * v[375];
	}
	else {
	};
	b2223 = b154;
	if (b2223) {
		v[2215] = v[2215] + v[2201] * v[540];
	}
	else {
	};
	v[2224] = v[2217];
	b2225 = b137;
	if (b2225) {
		v[2215] = v[2215] + v[142] * v[2224];
	}
	else {
	};
	v[2226] = v[2220];
	b2227 = b137;
	if (b2227) {
		v[2215] = v[2215] - v[142] * v[2226];
	}
	else {
	};
	v[4468] = v[2215] * v[444] + v[2202] * v[462] + v[200] * v[582];
	v[4469] = v[2215] * v[445] + v[2202] * v[463] + v[200] * v[583];
	v[4470] = v[2215] * v[446] + v[2202] * v[464] + v[200] * v[584];
	v[4471] = v[2215] * v[447] + v[2202] * v[465] + v[205] * v[582];
	v[4472] = v[2215] * v[448] + v[2202] * v[466] + v[205] * v[583];
	v[4473] = v[2215] * v[449] + v[2202] * v[467] + v[205] * v[584];
	v[4474] = v[2215] * v[450] + v[2202] * v[468] + v[211] * v[582];
	v[4475] = v[2215] * v[451] + v[2202] * v[469] + v[211] * v[583];
	v[4476] = v[2215] * v[452] + v[2202] * v[470] + v[211] * v[584];
	v[4477] = v[2215] * v[453] + v[2202] * v[471] - v[302] * v[582];
	v[4478] = v[2215] * v[454] + v[2202] * v[472] - v[302] * v[583];
	v[4479] = v[2215] * v[455] + v[2202] * v[473] - v[302] * v[584];
	v[4480] = v[2215] * v[456] + v[2202] * v[474] - v[307] * v[582];
	v[4481] = v[2215] * v[457] + v[2202] * v[475] - v[307] * v[583];
	v[4482] = v[2215] * v[458] + v[2202] * v[476] - v[307] * v[584];
	v[4483] = v[2215] * v[459] + v[2202] * v[477] - v[313] * v[582];
	v[4484] = v[2215] * v[460] + v[2202] * v[478] - v[313] * v[583];
	v[4485] = v[2215] * v[461] + v[2202] * v[479] - v[313] * v[584];
	for (i2194 = 1; i2194 <= 18; i2194++) {
		i2669 = (i2194 == 18 ? 1 : 0);
		i2668 = (i2194 == 17 ? 1 : 0);
		i2667 = (i2194 == 16 ? 1 : 0);
		i2666 = (i2194 == 15 ? 1 : 0);
		i2665 = (i2194 == 14 ? 1 : 0);
		i2664 = (i2194 == 13 ? 1 : 0);
		i2663 = (i2194 == 12 ? 1 : 0);
		i2662 = (i2194 == 11 ? 1 : 0);
		i2661 = (i2194 == 10 ? 1 : 0);
		i2660 = (i2194 == 9 ? 1 : 0);
		i2659 = (i2194 == 8 ? 1 : 0);
		i2658 = (i2194 == 7 ? 1 : 0);
		i2657 = (i2194 == 6 ? 1 : 0);
		i2656 = (i2194 == 5 ? 1 : 0);
		i2655 = (i2194 == 4 ? 1 : 0);
		i2654 = (i2194 == 3 ? 1 : 0);
		i2653 = (i2194 == 2 ? 1 : 0);
		i2652 = (i2194 == 1 ? 1 : 0);
		v[2287] = i2652 * v[444] + i2653 * v[445] + i2654 * v[446] + i2655 * v[447] + i2656 * v[448] + i2657 * v[449] + i2658 * v[450]
			+ i2659 * v[451] + i2660 * v[452] + i2661 * v[453] + i2662 * v[454] + i2663 * v[455] + i2664 * v[456] + i2665 * v[457]
			+ i2666 * v[458] + i2667 * v[459] + i2668 * v[460] + i2669 * v[461];
		v[2288] = i2652 * v[462] + i2653 * v[463] + i2654 * v[464] + i2655 * v[465] + i2656 * v[466] + i2657 * v[467] + i2658 * v[468]
			+ i2659 * v[469] + i2660 * v[470] + i2661 * v[471] + i2662 * v[472] + i2663 * v[473] + i2664 * v[474] + i2665 * v[475]
			+ i2666 * v[476] + i2667 * v[477] + i2668 * v[478] + i2669 * v[479];
		b2289 = b137;
		if (b2289) {
			v[2290] = -(v[142] * v[2287]);
		}
		else {
			v[2290] = 0e0;
		};
		v[2291] = v[2290];
		b2292 = b137;
		if (b2292) {
			v[2293] = v[142] * v[2287];
		}
		else {
			v[2293] = 0e0;
		};
		v[2294] = v[2293];
		v[2295] = 0e0;
		b2296 = b154;
		if (b2296) {
			v[2295] = v[2287] * v[540];
		}
		else {
		};
		v[2297] = 0e0;
		v[2298] = 0e0;
		b2299 = b157;
		if (b2299) {
			v[2295] = v[2295] + v[2287] * v[375];
			v[2298] = v[2201] * v[2287];
			v[2295] = v[2295] + v[2291] * v[2456];
			v[2297] = -(v[162] * v[2201] * v[2291]);
		}
		else {
		};
		v[2300] = 0e0;
		v[2301] = 0e0;
		b2302 = b157;
		if (b2302) {
			v[2301] = -(v[2287] * v[375]);
			v[2298] = -(v[2200] * v[2287]) + v[2298];
			v[2301] = v[2301] + v[2291] * v[2457];
			v[2300] = v[162] * v[2200] * v[2291];
			v[2291] = 0e0;
		}
		else {
		};
		v[2303] = 0e0;
		v[2304] = 0e0;
		b2305 = b164;
		if (b2305) {
			v[2301] = v[2301] + v[2287] * v[381];
			v[2304] = v[2200] * v[2287];
			v[2301] = v[2301] + v[2294] * v[2458];
			v[2303] = v[169] * v[2200] * v[2294];
		}
		else {
		};
		v[2306] = 0e0;
		b2307 = b164;
		if (b2307) {
			v[2306] = -(v[2287] * v[381]);
			v[2304] = -(v[2199] * v[2287]) + v[2304];
			v[2306] = v[2306] + v[2294] * v[2459];
			v[2297] = v[169] * v[2199] * v[2294] + v[2297];
			v[2294] = 0e0;
		}
		else {
		};
		b2308 = b171;
		if (b2308) {
			v[2306] = v[2306] + v[2287] * v[550];
			v[2287] = 0e0;
		}
		else {
		};
		b2309 = b239;
		if (b2309) {
			v[2310] = -(v[2288] * v[244]);
		}
		else {
			v[2310] = 0e0;
		};
		v[2311] = v[2310];
		b2312 = b239;
		if (b2312) {
			v[2313] = v[2288] * v[244];
		}
		else {
			v[2313] = 0e0;
		};
		v[2314] = v[2313];
		v[2315] = 0e0;
		b2316 = b256;
		if (b2316) {
			v[2315] = v[2288] * v[558];
		}
		else {
		};
		v[2317] = 0e0;
		v[2318] = 0e0;
		b2319 = b259;
		if (b2319) {
			v[2315] = v[2315] + v[2288] * v[394];
			v[2318] = v[2198] * v[2288];
			v[2315] = v[2315] + v[2311] * v[2462];
			v[2317] = -(v[2198] * v[2311] * v[264]);
		}
		else {
		};
		v[2320] = 0e0;
		v[2321] = 0e0;
		b2322 = b259;
		if (b2322) {
			v[2321] = -(v[2288] * v[394]);
			v[2318] = -(v[2197] * v[2288]) + v[2318];
			v[2321] = v[2321] + v[2311] * v[2463];
			v[2320] = v[2197] * v[2311] * v[264];
			v[2311] = 0e0;
		}
		else {
		};
		v[2323] = 0e0;
		v[2324] = 0e0;
		b2325 = b266;
		if (b2325) {
			v[2321] = v[2321] + v[2288] * v[400];
			v[2324] = v[2197] * v[2288];
			v[2321] = v[2321] + v[2314] * v[2464];
			v[2323] = v[2197] * v[2314] * v[271];
		}
		else {
		};
		v[2326] = 0e0;
		b2327 = b266;
		if (b2327) {
			v[2326] = -(v[2288] * v[400]);
			v[2324] = -(v[2196] * v[2288]) + v[2324];
			v[2326] = v[2326] + v[2314] * v[2465];
			v[2317] = v[2317] + v[2196] * v[2314] * v[271];
			v[2314] = 0e0;
		}
		else {
		};
		b2328 = b273;
		if (b2328) {
			v[2326] = v[2326] + v[2288] * v[568];
			v[2288] = 0e0;
		}
		else {
		};
		v[2359] = i2654 * v[200] + i2657 * v[205] + i2660 * v[211] + v[121] * v[2295] + v[124] * v[2301] + v[127] * v[2306]
			- v[223] * v[2315] - v[226] * v[2321] - v[229] * v[2326] - i2663 * v[302] - i2666 * v[307] - i2669 * v[313];
		v[2670] = v[13] * v[2359];
		v[2413] = v[2670] * v[360];
		v[2360] = i2653 * v[200] + i2656 * v[205] + i2659 * v[211] + v[120] * v[2295] + v[123] * v[2301] + v[126] * v[2306]
			- v[222] * v[2315] - v[225] * v[2321] - v[228] * v[2326] - i2662 * v[302] - i2665 * v[307] - i2668 * v[313];
		v[2671] = v[13] * v[2360];
		v[2382] = -(v[2671] * v[528]);
		v[2677] = v[2382] - v[2413] * v[359];
		v[2373] = v[2671] * v[359];
		v[2361] = i2652 * v[200] + i2655 * v[205] + i2658 * v[211] + v[119] * v[2295] + v[122] * v[2301] + v[125] * v[2306]
			- v[221] * v[2315] - v[224] * v[2321] - v[227] * v[2326] - i2661 * v[302] - i2664 * v[307] - i2667 * v[313];
		v[2672] = v[13] * v[2361];
		v[2383] = -(v[2672] * v[520]);
		v[2676] = v[2383] - v[2413] * v[358];
		v[2375] = v[2672] * v[358];
		v[2365] = v[2373] + v[2375];
		v[2675] = -(v[2365] * v[360]);
		v[2366] = v[13] * (v[2360] * v[323] + v[2361] * v[324]);
		v[2367] = v[13] * (v[2360] * v[326] + v[2361] * v[327]);
		v[2368] = v[13] * (v[2360] * v[329] + v[2361] * v[330]);
		v[2369] = v[13] * (v[2360] * v[332] + v[2361] * v[333]);
		v[2370] = v[13] * (v[2360] * v[335] + v[2361] * v[336]);
		v[2371] = v[13] * (v[2360] * v[338] + v[2361] * v[339]);
		v[2673] = v[200] * v[2366] + v[205] * v[2367] + v[211] * v[2368] - v[2369] * v[302] - v[2370] * v[307] - v[2371] * v[313];
		v[2372] = v[1008] * v[2365] + 2e0 * v[1008] * v[2413] + v[1024] * v[2670];
		v[2374] = v[13] * (v[1006] * v[2359] + v[1027] * v[2360]) + 2e0 * v[1023] * v[2373] + v[2673] * v[358];
		v[2376] = v[13] * (v[1005] * v[2359] + v[1027] * v[2361]) + 2e0 * v[1022] * v[2375] + v[2673] * v[359];
		v[2674] = ((v[2376] * v[341] + v[2374] * v[342] + v[2372] * v[343]) * v[356]) / v[351];
		v[2378] = v[2674] * v[343] + v[2372] * v[357];
		v[2678] = v[2378] + v[15] * (-v[2675] + v[2670] * v[536]);
		v[2379] = v[2674] * v[342] + v[2374] * v[357];
		v[2380] = v[2674] * v[341] + v[2376] * v[357];
		v[2381] = -(v[229] * v[2378]) - v[228] * v[2379] - v[227] * v[2380] + v[2371] * v[2479] + v[2383] * v[338]
			+ v[2382] * v[339] + v[2675] * v[340] - i2667 * v[582] - i2668 * v[583] - i2669 * v[584] + v[2670] * v[979];
		v[2384] = -(v[226] * v[2378]) - v[225] * v[2379] - v[224] * v[2380] + v[2370] * v[2479] + v[2383] * v[335]
			+ v[2382] * v[336] + v[2675] * v[337] - i2664 * v[582] - i2665 * v[583] - i2666 * v[584] + v[2670] * v[960];
		v[2385] = -(v[223] * v[2378]) - v[222] * v[2379] - v[221] * v[2380] + v[2369] * v[2479] + v[2383] * v[332]
			+ v[2382] * v[333] + v[2675] * v[334] - i2661 * v[582] - i2662 * v[583] - i2663 * v[584] + v[2670] * v[941];
		v[2386] = v[127] * v[2378] + v[126] * v[2379] + v[125] * v[2380] - v[2368] * v[2479] - v[2383] * v[329] - v[2382] * v[330]
			- v[2675] * v[331] + i2658 * v[582] + i2659 * v[583] + i2660 * v[584] + v[2670] * v[922];
		v[2387] = v[124] * v[2378] + v[123] * v[2379] + v[122] * v[2380] - v[2367] * v[2479] - v[2383] * v[326] - v[2382] * v[327]
			- v[2675] * v[328] + i2655 * v[582] + i2656 * v[583] + i2657 * v[584] + v[2670] * v[900];
		v[2388] = v[121] * v[2378] + v[120] * v[2379] + v[119] * v[2380] - v[2366] * v[2479] - v[2383] * v[323] - v[2382] * v[324]
			- v[2675] * v[325] + i2652 * v[582] + i2653 * v[583] + i2654 * v[584] + v[2670] * v[869];
		b2389 = b273;
		if (b2389) {
			v[2317] = v[2317] + v[2381] * v[568];
		}
		else {
		};
		v[2390] = 0e0;
		b2391 = b266;
		if (b2391) {
			v[2390] = v[2381] * v[2465];
			v[2317] = v[2317] - v[2381] * v[400];
		}
		else {
		};
		b2392 = b266;
		if (b2392) {
			v[2323] = v[2323] - v[2384] * v[400];
			v[2317] = v[2317] - v[2323];
			v[2390] = v[2390] + v[271] * (-v[2324] + v[2384] * v[397]);
		}
		else {
		};
		v[2393] = 0e0;
		b2394 = b259;
		if (b2394) {
			v[2320] = v[2320] - v[2384] * v[394];
			v[2393] = v[2384] * v[2463];
			v[2317] = v[2317] + v[2320];
		}
		else {
		};
		b2395 = b259;
		if (b2395) {
			v[2317] = v[2317] + v[2385] * v[394];
			v[2393] = v[2393] + v[264] * (-v[2318] + v[2385] * v[730]);
		}
		else {
		};
		b2396 = b256;
		if (b2396) {
			v[2317] = v[2317] + v[2385] * v[558];
		}
		else {
		};
		v[2397] = v[2390];
		b2398 = b239;
		if (b2398) {
			v[2317] = v[2317] + v[2397] * v[244];
		}
		else {
		};
		v[2399] = v[2393];
		b2400 = b239;
		if (b2400) {
			v[2317] = v[2317] - v[2399] * v[244];
		}
		else {
		};
		b2401 = b171;
		if (b2401) {
			v[2297] = v[2297] + v[2386] * v[550];
		}
		else {
		};
		v[2402] = 0e0;
		b2403 = b164;
		if (b2403) {
			v[2402] = v[2386] * v[2459];
			v[2297] = v[2297] - v[2386] * v[381];
		}
		else {
		};
		b2404 = b164;
		if (b2404) {
			v[2303] = v[2303] - v[2387] * v[381];
			v[2297] = v[2297] - v[2303];
			v[2402] = v[2402] + v[169] * (-v[2304] + v[2387] * v[378]);
		}
		else {
		};
		v[2405] = 0e0;
		b2406 = b157;
		if (b2406) {
			v[2300] = v[2300] - v[2387] * v[375];
			v[2405] = v[2387] * v[2457];
			v[2297] = v[2297] + v[2300];
		}
		else {
		};
		b2407 = b157;
		if (b2407) {
			v[2297] = v[2297] + v[2388] * v[375];
			v[2405] = v[2405] + v[162] * (-v[2298] + v[2388] * v[745]);
		}
		else {
		};
		b2408 = b154;
		if (b2408) {
			v[2297] = v[2297] + v[2388] * v[540];
		}
		else {
		};
		v[2409] = v[2402];
		b2410 = b137;
		if (b2410) {
			v[2297] = v[2297] + v[142] * v[2409];
		}
		else {
		};
		v[2411] = v[2405];
		b2412 = b137;
		if (b2412) {
			v[2297] = v[2297] - v[142] * v[2411];
		}
		else {
		};
		v[4490] = v[200] * v[2380] + v[2297] * v[444] + v[2317] * v[462] + v[15] * (-(v[200] * v[2676]) + v[2671] * v[522])
			+ v[2295] * v[582];
		v[4491] = v[200] * v[2379] + v[2297] * v[445] + v[2317] * v[463] + v[15] * (-(v[200] * v[2677]) + v[2672] * v[522])
			+ v[2295] * v[583];
		v[4492] = v[200] * v[2678] + v[2297] * v[446] + v[2317] * v[464] + v[2295] * v[584];
		v[4493] = v[205] * v[2380] + v[2297] * v[447] + v[2317] * v[465] + v[15] * (-(v[205] * v[2676]) + v[2671] * v[523])
			+ v[2301] * v[582];
		v[4494] = v[205] * v[2379] + v[2297] * v[448] + v[2317] * v[466] + v[15] * (-(v[205] * v[2677]) + v[2672] * v[523])
			+ v[2301] * v[583];
		v[4495] = v[205] * v[2678] + v[2297] * v[449] + v[2317] * v[467] + v[2301] * v[584];
		v[4496] = v[211] * v[2380] + v[2297] * v[450] + v[2317] * v[468] + v[15] * (-(v[211] * v[2676]) + v[2671] * v[524])
			+ v[2306] * v[582];
		v[4497] = v[211] * v[2379] + v[2297] * v[451] + v[2317] * v[469] + v[15] * (-(v[211] * v[2677]) + v[2672] * v[524])
			+ v[2306] * v[583];
		v[4498] = v[211] * v[2678] + v[2297] * v[452] + v[2317] * v[470] + v[2306] * v[584];
		v[4499] = -(v[2380] * v[302]) + v[2297] * v[453] + v[2317] * v[471] + v[15] * (v[2676] * v[302] + v[2671] * v[525])
			- v[2315] * v[582];
		v[4500] = -(v[2379] * v[302]) + v[2297] * v[454] + v[2317] * v[472] + v[15] * (v[2677] * v[302] + v[2672] * v[525])
			- v[2315] * v[583];
		v[4501] = -(v[2678] * v[302]) + v[2297] * v[455] + v[2317] * v[473] - v[2315] * v[584];
		v[4502] = -(v[2380] * v[307]) + v[2297] * v[456] + v[2317] * v[474] + v[15] * (v[2676] * v[307] + v[2671] * v[526])
			- v[2321] * v[582];
		v[4503] = -(v[2379] * v[307]) + v[2297] * v[457] + v[2317] * v[475] + v[15] * (v[2677] * v[307] + v[2672] * v[526])
			- v[2321] * v[583];
		v[4504] = -(v[2678] * v[307]) + v[2297] * v[458] + v[2317] * v[476] - v[2321] * v[584];
		v[4505] = -(v[2380] * v[313]) + v[2297] * v[459] + v[2317] * v[477] + v[15] * (v[2676] * v[313] + v[2671] * v[527])
			- v[2326] * v[582];
		v[4506] = -(v[2379] * v[313]) + v[2297] * v[460] + v[2317] * v[478] + v[15] * (v[2677] * v[313] + v[2672] * v[527])
			- v[2326] * v[583];
		v[4507] = -(v[2678] * v[313]) + v[2297] * v[461] + v[2317] * v[479] - v[2326] * v[584];
		Rc[i2194 - 1] += v[4467 + i2194];
		for (i2267 = 1; i2267 <= 18; i2267++) {
			Kc[i2194 - 1][i2267 - 1] += v[4489 + i2267];
		};/* end for */
	};/* end for */
	//fn[0] = (*epsn)*v[212];
	//fn[1] = (*epsn)*v[213];
	//fn[2] = (*epsn)*v[214];
	//PrintPtr(fn, 3);
	//PrintPtr(Rc, 18);


#pragma endregion
}