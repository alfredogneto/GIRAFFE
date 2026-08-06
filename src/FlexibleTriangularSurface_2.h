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
	
	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calculates diversos vectors associated the surface for a pair of convective coordinates
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Data the point xS, calculates the coordinates (zeta,theta) referring to the minimum distance
	void FillNodes();														//Updates the variables internal of the surface, for grab info of the pilot node for use subsequent with position updated
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calculates contributions of contact between sphere and surface
	void ContactSphereSurfaceSliding (double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calculates contributions of contact between sphere and surface
	void SaveConfiguration();												//Saves vectors of configuration converged
	void CenterPoint(Matrix* center);										//Returns coordinates global of the point central of the surface the be used for calculations coarse of its location (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Performs initial guess for the variables zeta and theta
	bool Check();															//Checks inconsistencies for avoid errors of execution
	void NormalExt(double* zeta, double* theta, Matrix* n);							//Normal exterior the surface in the position chosen
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtains point of the surface
	void SetMinMaxRange();

	void UpdateBox();				//Updates bounding box
	//Variables internal
	double offset;		//offset in the normal direction the be applied in the surface, in relação the equação of the parametrização proposta

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

