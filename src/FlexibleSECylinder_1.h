#pragma once
#include "Surface.h"

class FlexibleSECylinder_1 :
	public Surface
{
public:
	FlexibleSECylinder_1();
	~FlexibleSECylinder_1();
	void PreCalc();
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calcula diversos vetores associados à superfície para um par de coordenadas convectivas
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes à mínima distância
	void FillNodes();														//Atualiza as variáveis internas da superfície, para pegarem info do pilot node para uso posterior com posição atualizada
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribuições de contato entre esfera e superfície
	void ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribuições de contato entre esfera e superfície
	void SaveConfiguration();												//Salva vetores de configuração convergida
	void CenterPoint(Matrix* center);										//Retorna coordenadas globais do ponto central da superfície a ser utilizado para cálculos grosseiros de sua localização (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Realiza chute inicial para as variáveis zeta e theta
	bool Check();															//Checa inconsistências para evitar erros de execução
	void Gamma(double* zeta, double* theta, Matrix* G);						//Posição da superfície nas coordenadas convectivas escolhidas
	void NormalExt(double* zeta, double* theta, Matrix* n);					//Normal exterior à superfície na posição escolhida
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtem ponto da superficie
	void SetMinMaxRange();

	void GammaIntersectionLine(Matrix* p_PO, Matrix* p_t, double* p_lambda, Matrix* p_c);	//Intersecção Gamma com a reta
	double Curvature(double zeta, double theta, Matrix* e_dir);			//Retorna a curvatura da superfície segundo certa direção e calculada em certo conjunto de coordenadas convectivas
	void UpdateBox();				//Atualiza bounding box
	//Variáveis internas
	double a, b, n, e;
	int csA, csB;
	bool flag_normal_int;
	Matrix* vNR;

	Matrix* x_AAi;											
	Matrix* x_BAi;											
	Matrix* Q_AAi;													
	Matrix* Q_BAi;
	Matrix* Q_AAic;
	Matrix* Q_BAic;
	Matrix* d_A;
	Matrix* dui_A;
	Matrix* ddui_A;
	Matrix* alpha_AA;
	Matrix* alpha_BA;
	Matrix* Q0A;
	Matrix* Q0B;

	Matrix* x_AAp;
	Matrix* x_BAp;
	Matrix* Q_AAp;
	Matrix* Q_BAp;

	double** aQ_AAi;
	double** aQ_BAi;

	Matrix* I3;																//Identidade de ordem 3
};

