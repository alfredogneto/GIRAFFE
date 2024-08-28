#pragma once
#include "Surface.h"

class RigidOscillatorySurface_1 :
	public Surface
{
public:
	RigidOscillatorySurface_1();
	~RigidOscillatorySurface_1();
	void PreCalc();
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f);
	
	int* points;

	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calcula diversos vetores associados � superf�cie para um par de coordenadas convectivas
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes � m�nima dist�ncia
	void FillNodes();														//Atualiza as vari�veis internas da superf�cie, para pegarem info do pilot node para uso posterior com posi��o atualizada
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribui��es de contato entre esfera e superf�cie
	void ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribui��es de contato entre esfera e superf�cie
	void SaveConfiguration();												//Salva vetores de configura��o convergida
	void CenterPoint(Matrix* center);										//Retorna coordenadas globais do ponto central da superf�cie a ser utilizado para c�lculos grosseiros de sua localiza��o (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Realiza chute inicial para as vari�veis zeta e theta
	void Gamma(Matrix *Gamma, double *zeta, double *theta);					//Retorna em G o valor de Gamma para zeta e theta (tomando posi��o atual)
	bool Check();															//Checa inconsist�ncias para evitar erros de execu��o
	void NormalExt(double* zeta, double* theta, Matrix* n);							//Normal exterior � superf�cie na posi��o escolhida
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtem ponto da superficie
	void UpdateBox();				//Atualiza bounding box
	void SetMinMaxRange();

	//Vari�veis internas
	Matrix* xP_i;
	Matrix* e1_i;
	Matrix* e2_i;

	Matrix* xP_p;
	Matrix* e1_p;
	Matrix* e2_p;

	int number_CS;
	double A_1, A_2, A_12, lambda_1, lambda_2, phi_1, phi_2;
	double waves_1, waves_2;


	Matrix* I3;
	Matrix* vNR;
};
