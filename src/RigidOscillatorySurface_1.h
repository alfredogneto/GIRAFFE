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

	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calculates diversos vectors associated the surface for a pair of convective coordinates
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Data the point xS, calculates the coordinates (zeta,theta) referring to the minimum distance
	void FillNodes();														//Updates the variables internal of the surface, for grab info of the pilot node for use subsequent with position updated
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calculates contributions of contact between sphere and surface
	void ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calculates contributions of contact between sphere and surface
	void SaveConfiguration();												//Saves vectors of configuration converged
	void CenterPoint(Matrix* center);										//Returns coordinates global of the point central of the surface the be used for calculations coarse of its location (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Performs initial guess for the variables zeta and theta
	void Gamma(Matrix *Gamma, double *zeta, double *theta);					//Returns in G the value of Gamma for zeta and theta (tomando position current)
	bool Check();															//Checks inconsistencies for avoid errors of execution
	void NormalExt(double* zeta, double* theta, Matrix* n);							//Normal exterior the surface in the position chosen
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtains point of the surface
	void UpdateBox();				//Updates bounding box
	void SetMinMaxRange();

	//Variables internal
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
