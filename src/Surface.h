#pragma once
#include <stdio.h>
#include "Box.h"

class Matrix;
class NSContactData;

class Surface
{
public:
	Surface();
	virtual ~Surface();

	int *nodes;		//Nós globais - conectividade
	int number;		//ID da superficie
	int n_nodes;	//Numero de nós da superficie
	int **DOFs;		//Indica para a indexação de cada grau de liberdade, 1 ou 0, ativo ou inativo para o elemento em questão
	int nDOFs;		//Numero de GL
	int VTK_type;	//Tipo de celula para VTK
	int *VTK_nodes;	//Indexação para converter numeração do formato giraffe para o formato da celula equivalente do paraview
	int **GLs;		//Ponteiro para os GL globais utilizados na superficie
	int pilot_node;											//Nó piloto usado para a construção da parametrização
	bool pilot_is_used;										//flag que indica se o pilot node e usado

	Box box;		//Bounding box

	//Degeneration (true or false for first and second convective coordinates)
	bool degeneration[2];
	//Degenerated coordinates (only considered if boolean degeneration values are TRUE)
	double ***deg_coordinates;
	//Divisions on degeneration
	int div1, div2;
	//Bool indicators of coordinates
	bool entered_u1, entered_u2;
	double u1_min, u1_max, u2_min, u2_max, u1_range, u2_range;
	bool alloced_degeneration;
	int alloced_div1, alloced_div2;

	virtual void PreCalc() = 0;
	virtual bool Read(FILE *f) = 0;
	virtual void SetMinMaxRange() = 0;

	bool ReadCommon(FILE *f);
	void AllocDegeneration();
	void FreeDegeneration();
	void DegenerationPreCalc();
	void InitializeDegeneration();

	virtual void Write(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual bool Check() = 0;															//Checa inconsistências para evitar erros de execução
	//Funções virtuais para utilização do contato do tipo NSSS
	virtual void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp) = 0;				//Calcula diversos vetores associados a superficie para um par de coordenadas convectivas
	virtual void FindMinimimumParameters(Matrix* xS, NSContactData* cd) = 0;			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes a minima distancia
	virtual void FillNodes() = 0;														//Atualiza as variaveis internas da superficie, para pegarem info do pilot node para uso posterior com posição atualizada
	virtual void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius) = 0;		//Calcula contribuições de contato entre esfera e superficie
	virtual void ContactSphereSurfaceSliding (double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius) = 0;		//Calcula contribuições de contato entre esfera e superficie
	virtual void SaveConfiguration() = 0;												//Salva vetores de configuração convergida
	virtual void CenterPoint(Matrix* center) = 0;										//Retorna coordenadas globais do ponto central da superficie a ser utilizado para calculos grosseiros de sua localização (pinball)
	virtual void InitialGuess(Matrix* xS, double** convective, int n_solutions) = 0;	//Realiza chute inicial para as variaveis zeta e theta
	virtual void UpdateBox() = 0;														//Obtem bounding box
	virtual void NormalExt(double* zeta, double* theta, Matrix* n) = 0;					//Obtem a normal exterior da superficie
	virtual void SurfacePoint(double& zeta, double& theta, Matrix& point) = 0;			//Obtem ponto da superficie
};

