#pragma once
#include <stdio.h>
#include "Surface.h"

#include "Matrix.h"

class NURBSSurface;

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

	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp);				//Calculates diversos vectors associated the surface for a pair of convective coordinates
	void FindMinimimumParameters(Matrix* xS, NSContactData* cd);			//Data the point xS, calculates the coordinates (zeta,theta) referring to the minimum distance
	void FillNodes();														//Updates the variables internal of the surface, for grab info of the pilot node for use subsequent with position updated
	void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calculates contributions of contact between sphere and surface
	void ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius);		//Calculates contributions of contact between sphere and surface
	void SaveConfiguration();												//Saves vectors of configuration converged
	void CenterPoint(Matrix* center);										//Returns coordinates global of the point central of the surface the be used for calculations coarse of its location (pinball)
	void InitialGuess(Matrix* xS, double** convective, int n_solutions);	//Performs initial guess for the variables zeta and theta
	bool Check();															//Checks inconsistencies for avoid errors of execution
	void NormalExt(double* zeta, double* theta, Matrix* n);					//Normal exterior the surface in the position chosen
	void SurfacePoint(double& zeta, double& theta, Matrix& point);			//Obtains point of the surface
	void SetMinMaxRange();

	void UpdateBox();														//Updates bounding box
	
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

