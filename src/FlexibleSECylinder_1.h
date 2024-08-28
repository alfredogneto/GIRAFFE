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

	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calcula diversos vetores associados � superf�cie para um par de coordenadas convectivas
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes � m�nima dist�ncia
	void FillNodes();														//Atualiza as vari�veis internas da superf�cie, para pegarem info do pilot node para uso posterior com posi��o atualizada
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribui��es de contato entre esfera e superf�cie
	void ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribui��es de contato entre esfera e superf�cie
	void SaveConfiguration();												//Salva vetores de configura��o convergida
	void CenterPoint(Matrix* center);										//Retorna coordenadas globais do ponto central da superf�cie a ser utilizado para c�lculos grosseiros de sua localiza��o (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Realiza chute inicial para as vari�veis zeta e theta
	bool Check();															//Checa inconsist�ncias para evitar erros de execu��o
	void Gamma(double* zeta, double* theta, Matrix* G);						//Posi��o da superf�cie nas coordenadas convectivas escolhidas
	void NormalExt(double* zeta, double* theta, Matrix* n);					//Normal exterior � superf�cie na posi��o escolhida
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtem ponto da superficie
	void SetMinMaxRange();

	void GammaIntersectionLine(Matrix* p_PO, Matrix* p_t, double* p_lambda, Matrix* p_c);	//Intersec��o Gamma com a reta
	double Curvature(double zeta, double theta, Matrix* e_dir);			//Retorna a curvatura da superf�cie segundo certa dire��o e calculada em certo conjunto de coordenadas convectivas
	void UpdateBox();				//Atualiza bounding box
	//Vari�veis internas
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

