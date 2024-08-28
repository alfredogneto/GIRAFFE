#pragma once
#include "Surface.h"

class FlexibleTriangularSurface_2 :
	public Surface
{
public:
	FlexibleTriangularSurface_2();
	~FlexibleTriangularSurface_2();
	void PreCalc();
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f);
	
	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calcula diversos vetores associados � superf�cie para um par de coordenadas convectivas
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes � m�nima dist�ncia
	void FillNodes();														//Atualiza as vari�veis internas da superf�cie, para pegarem info do pilot node para uso posterior com posi��o atualizada
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribui��es de contato entre esfera e superf�cie
	void ContactSphereSurfaceSliding (double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribui��es de contato entre esfera e superf�cie
	void SaveConfiguration();												//Salva vetores de configura��o convergida
	void CenterPoint(Matrix* center);										//Retorna coordenadas globais do ponto central da superf�cie a ser utilizado para c�lculos grosseiros de sua localiza��o (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Realiza chute inicial para as vari�veis zeta e theta
	bool Check();															//Checa inconsist�ncias para evitar erros de execu��o
	void NormalExt(double* zeta, double* theta, Matrix* n);							//Normal exterior � superf�cie na posi��o escolhida
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtem ponto da superficie
	void SetMinMaxRange();

	void UpdateBox();				//Atualiza bounding box
	//Vari�veis internas
	double offset;		//offset na dire��o normal a ser aplicado na superf�cie, em rela��o � equa��o da parametriza��o proposta

	Matrix* xA_p;
	Matrix* xB_p;
	Matrix* xC_p;
	Matrix* xD_p;
	Matrix* xE_p;
	Matrix* xF_p;

	Matrix* xA_i;
	Matrix* xB_i;
	Matrix* xC_i;
	Matrix* xD_i;
	Matrix* xE_i;
	Matrix* xF_i;

	Matrix* d_A;
	Matrix* dui_A;
	Matrix* ddui_A;

	Matrix* vNR;
};

