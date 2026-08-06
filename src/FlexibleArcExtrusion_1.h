#pragma once
#include "Surface.h"

class FlexibleArcExtrusion_1 :
	public Surface
{
public:
	FlexibleArcExtrusion_1();
	~FlexibleArcExtrusion_1();
	void PreCalc();
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zetai, double* thi, double* zetap, double* thp);				//Calculates diversos vectors associated the surface for a pair of convective coordinates
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
	void UpdateBox();														//Obtains bounding box
	void SetMinMaxRange();

	//Variables of the Element
	int arc_ID;													//ID of the arc the be extrudado
	int cs;                                                     //ID of the coordinate system for posicionar the arc in the espaço
	double* radius;                                             //Radius of curvature of the arc
	Matrix* c_point;											//Center of curvature of the arc
	Matrix* i_point;											//Point initial of the arc
	Matrix* f_point;											//Point final of the arc

	double* theta_i;
	double* theta_f;

	//Variables for escolha of the lado concavo or convexo
	bool flag_normal_int;									//Flag for indicate use of normal interior/exterior
	
	Matrix* x_AAi;
	Matrix* x_BAi;
	Matrix* Q_AAi;
	Matrix* Q_BAi;
	Matrix* Q_AAic;
	Matrix* Q_BAic;
	Matrix* d_A;
	Matrix* dui_A;
	Matrix* ddui_A;
	Matrix* alpha_AA;
	Matrix* alpha_BA;
	Matrix* Q0A;

	Matrix* x_AAp;
	Matrix* x_BAp;
	Matrix* Q_AAp;
	Matrix* Q_BAp;

	double** aQ_AAi;
	double** aQ_BAi;

	Matrix* I3;																//Identidade of order 3
};

