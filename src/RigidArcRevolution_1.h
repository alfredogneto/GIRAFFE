#pragma once
#include "Surface.h"

class RigidArcRevolution_1 :
	public Surface
{
public:
	RigidArcRevolution_1();
	~RigidArcRevolution_1();
	void PreCalc();
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* phii, double* thi, double* phip, double* thp);				//Calcula diversos vetores associados a superficie para um par de coordenadas convectivas
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes a minima distancia
	void FillNodes();														//Atualiza as variaveis internas da superficie, para pegarem info do pilot node para uso posterior com posição atualizada
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribuições de contato entre esfera e superficie
	void ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calcula contribuições de contato entre esfera e superficie
	void SaveConfiguration();												//Salva vetores de configuração convergida
	void CenterPoint(Matrix* center);										//Retorna coordenadas globais do ponto central da superficie a ser utilizado para calculos grosseiros de sua localização (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Realiza chute inicial para as variaveis zeta e theta
	bool Check();															//Checa inconsistências para evitar erros de execução
	void NormalExt(double* zeta, double* theta, Matrix* n);							//Normal exterior a superficie na posição escolhida
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtem ponto da superficie
	void SetMinMaxRange();

	void UpdateBox();														//Obtem bounding box


	//Variaveis do Elemento
	int arc_ID;													//ID do arco a ser extrudado
	int cs;                                                     //ID do sistema de coordenadas para posicionar o arco no espaço
	double rev_Ang;                                             //Angle of revolution - default = 2*pi
	double x_fac;												//Fator para representar ovalização no semi-eixo x - default = 1
	double z_fac;												//Fator para representar ovalização no semi-eixo z - default = 1

	double* radius;                                              //Raio de curvatura do arco
	Matrix* c_point;											//Centro de curvatura do arco
	Matrix* i_point;											//Ponto inicial do arco
	Matrix* f_point;											//Ponto final do arco

	double* theta_i;
	double* theta_f;

	//Variaveis para escolha do lado concavo ou convexo
	bool flag_normal_int;									//Flag para indicar uso de normal interior/exterior


	Matrix* x_ABi;
	Matrix* Q_ABi;
	Matrix* Q_ABic;
	Matrix* d_B;
	Matrix* dui_B;
	Matrix* ddui_B;
	Matrix* alpha_AB;
	Matrix* Q0B;

	Matrix* x_ABp;
	Matrix* Q_ABp;

	double** aQ_ABi;
	
	Matrix* I3;																//Identidade de ordem 3
};

