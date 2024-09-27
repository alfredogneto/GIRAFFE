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
	
	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calcula diversos vetores associados a superficie para um par de coordenadas convectivas
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes a minima distancia
	void FillNodes();														//Atualiza as variaveis internas da superficie, para pegarem info do pilot node para uso posterior com posição atualizada
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribuições de contato entre esfera e superficie
	void ContactSphereSurfaceSliding (double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribuições de contato entre esfera e superficie
	void SaveConfiguration();												//Salva vetores de configuração convergida
	void CenterPoint(Matrix* center);										//Retorna coordenadas globais do ponto central da superficie a ser utilizado para calculos grosseiros de sua localização (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Realiza chute inicial para as variaveis zeta e theta
	bool Check();															//Checa inconsistências para evitar erros de execução
	void NormalExt(double* zeta, double* theta, Matrix* n);							//Normal exterior a superficie na posição escolhida
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtem ponto da superficie
	void SetMinMaxRange();

	void UpdateBox();				//Atualiza bounding box
	//Variaveis internas
	double offset;		//offset na direção normal a ser aplicado na superficie, em relação a equação da parametrização proposta

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

