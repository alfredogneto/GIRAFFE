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

	int *nodes;		//Nodes global - conectividade
	int number;		//ID of the surface
	int n_nodes;	//Number of nodes of the surface
	int **DOFs;		//Indicates for the indexing of each degree of freedom, 1 or 0, active or inactive for the element in question
	int nDOFs;		//Number of GL
	int VTK_type;	//Type of celula for VTK
	int *VTK_nodes;	//Indexing for convert numbering of the formato giraffe for the formato of the celula equivalente of the paraview
	int **GLs;		//Pointer for the GL global used in the surface
	int pilot_node;											//Node piloto used for the construction of the parametrização
	bool pilot_is_used;										//flag that indicates if the pilot node and used

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
	virtual bool Check() = 0;															//Checks inconsistencies for avoid errors of execution
	//Functions virtuais for utilização of the contact of the type NSSS
	virtual void Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp) = 0;				//Calculates diversos vectors associated the surface for a pair of convective coordinates
	virtual void FindMinimimumParameters(Matrix* xS, NSContactData* cd) = 0;			//Data the point xS, calculates the coordinates (zeta,theta) referring to the minimum distance
	virtual void FillNodes() = 0;														//Updates the variables internal of the surface, for grab info of the pilot node for use subsequent with position updated
	virtual void ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius) = 0;		//Calculates contributions of contact between sphere and surface
	virtual void ContactSphereSurfaceSliding (double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius) = 0;		//Calculates contributions of contact between sphere and surface
	virtual void SaveConfiguration() = 0;												//Saves vectors of configuration converged
	virtual void CenterPoint(Matrix* center) = 0;										//Returns coordinates global of the point central of the surface the be used for calculations coarse of its location (pinball)
	virtual void InitialGuess(Matrix* xS, double** convective, int n_solutions) = 0;	//Performs initial guess for the variables zeta and theta
	virtual void UpdateBox() = 0;														//Obtains bounding box
	virtual void NormalExt(double* zeta, double* theta, Matrix* n) = 0;					//Obtains the normal exterior of the surface
	virtual void SurfacePoint(double& zeta, double& theta, Matrix& point) = 0;			//Obtains point of the surface
};

