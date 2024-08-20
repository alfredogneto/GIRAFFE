#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Surface.h"
#include "NURBSSurface.h"

class RigidNURBS_1 :
	public Surface
{
public:
	RigidNURBS_1();
	~RigidNURBS_1();

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
	void NormalExt(double* zeta, double* theta, Matrix* n);					//Normal exterior à superfície na posição escolhida
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtem ponto da superficie
	void SetMinMaxRange();

	void UpdateBox();														//Atualiza bounding box
	
	//Internal and storing variables
	int number_CS;
	int CADData_ID;

	bool invert_normal;
	Matrix Qi;						//Rotation "i"
	Matrix Qip;						//Rotation "i+1"
	Matrix xi;						//Pilot position "i"
	Matrix xip;						//Pilot position "i+1"
	Matrix I3;

	//Pointer to NURBS Surface
	NURBSSurface* surf;
};

