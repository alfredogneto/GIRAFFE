#pragma once
#include "Particle.h"

class Matrix;
class MatrixFloat;

class VEMPolyhedron :
	public Particle
{
public:
	VEMPolyhedron();
	virtual ~VEMPolyhedron();

	int CADDATA_ID;								//CAD ID associated
	
	//float bv_factor;		//Controls the size of bounding volumes of edges and vertices
	float inc_len_factor;	//Controls inflation of bounding volumes

	//Variables for kinematics
	Matrix** vertices_i;
	Matrix** vertices_p;
	int n_vertices;
	int n_faces;
	int n_edges;
	Matrix* x0i;
	Matrix* x0p;

	//To coumpound the convergence criteria
	/*int*		count_local_abs_load;
	double*		local_abs_load;*/

	double*		local_residual;		//residual to compound the discretized weak form
	double*		local_e_load;		//residual due to external loads (volume)
	double**	local_stiffness;	//stiffness to compound the discretized linearization of the weak form
	double**	local_mass;			//mass matrix
	double**	local_damping;		//Damping matrix

	//Explicit
	void EvaluateExplicit();
	void EvaluateAccelerations();
	void InitialEvaluations();

	//Geometric entities
	float deformed_volume;
	MatrixFloat* deformed_barycenter;
	float deformed_radius;
	float* deformed_faces_radii;
	float* deformed_edges_lengths;
	void UpdateGeometricEntities();
	void UpdateVolumeAndBarycenter();

	double mass;					//particle total mass
	//Load factor
	double l_factor;
	
	bool allocced;
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs);
	bool Check();
	void Alloc();
	void Free();
	void Mount();
	void MountGlobal();
	void MountFieldLoading();
	void InertialContributions();
	void DampingContributions();
	void UpdateVariables();
	void PreCalc();								//Pré-cálculo de variáveis que é feito uma única vez no início
	void UpdateBoundingVolumes();
	void SaveLagrange();						//Salva variáveis
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);
	void PostProcessing();
	
	
	bool damping_flag;				//Flag to turn damping on/off
	bool first_damping_evaluation;	//Flag first damping
	double damping_factor;			//Desired damping ratio for the first eigenvalue
	double alpha, beta;				//Rayleigh damping factors

	//AceGen Variables (pointers - when possible)
	double *Xref;															//Reference coordinates
	double *Xcur;															//Current coordinates
	double *displacement;													//Displacements
	double *domainData;														//Domain data
	double *elementData;													//Element data
	double *elementPostprocessing;											//Post processing data
	int* triangles;															//Connectivity for AceGen
	int* tetrahedrons;														//Tetrahedrons - only used when choosing stabilization type 1
	double* ht;																//History variable
	double* hp;																//History variable
	double *residual;														//Residual
	double *residualload;													//Residual for load (volumetric)
	double**tangent;														//Tangent
	double**massmatrix;														//mass matrix
	//Stabilization variables
	int stab_flag;
	double stab_factor;
	double stab_factor_mass;
	double stab_factor_load;
	bool lumped_mass_flag;
};

