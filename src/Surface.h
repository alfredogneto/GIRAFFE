#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "NSContactData.h"
#include "Box.h"

class Surface
{
public:
	Surface();
	virtual ~Surface();

	int *nodes;		//N�s globais - conectividade
	int number;		//ID da superf�cie
	int n_nodes;	//N�mero de n�s da superf�cie
	int **DOFs;		//Indica para a indexa��o de cada grau de liberdade, 1 ou 0, ativo ou inativo para o elemento em quest�o
	int nDOFs;		//N�mero de GL
	int VTK_type;	//Tipo de c�lula para VTK
	int *VTK_nodes;	//Indexa��o para converter numera��o do formato giraffe para o formato da c�lula equivalente do paraview
	int **GLs;		//Ponteiro para os GL globais utilizados na superf�cie
	int pilot_node;											//N� piloto usado para a constru��o da parametriza��o
	bool pilot_is_used;										//flag que indica se o pilot node � usado

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
	virtual bool Check() = 0;															//Checa inconsist�ncias para evitar erros de execu��o
	//Fun��es virtuais para utiliza��o do contato do tipo NSSS
	virtual void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp) = 0;				//Calcula diversos vetores associados � superf�cie para um par de coordenadas convectivas
	virtual void FindMinimimumParameters(Matrix* xS, NSContactData* cd) = 0;			//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes � m�nima dist�ncia
	virtual void FillNodes() = 0;														//Atualiza as vari�veis internas da superf�cie, para pegarem info do pilot node para uso posterior com posi��o atualizada
	virtual void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius) = 0;		//Calcula contribui��es de contato entre esfera e superf�cie
	virtual void ContactSphereSurfaceSliding (double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius) = 0;		//Calcula contribui��es de contato entre esfera e superf�cie
	virtual void SaveConfiguration() = 0;												//Salva vetores de configura��o convergida
	virtual void CenterPoint(Matrix* center) = 0;										//Retorna coordenadas globais do ponto central da superf�cie a ser utilizado para c�lculos grosseiros de sua localiza��o (pinball)
	virtual void InitialGuess(Matrix* xS, double** convective, int n_solutions) = 0;	//Realiza chute inicial para as vari�veis zeta e theta
	virtual void UpdateBox() = 0;														//Obtem bounding box
	virtual void NormalExt(double* zeta, double* theta, Matrix* n) = 0;					//Obtem a normal exterior da superficie
	virtual void SurfacePoint(double& zeta, double& theta, Matrix& point) = 0;			//Obtem ponto da superficie
};
