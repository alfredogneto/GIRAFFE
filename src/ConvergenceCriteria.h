#pragma once
#include <stdio.h>

class ConvergenceCriteria
{
public:
	ConvergenceCriteria();
	~ConvergenceCriteria();

	//Functions for estabelecer and check convergence
	void EstablishResidualCriteria();
	bool CheckResidualConvergence();
	bool CheckGLConvergence();

	bool Read(FILE *f);			//Reading
	void Write(FILE *f);		//Recording

	double force_tol;			//Factor applied the norm of the force vector of the system for estabelecer convergence criterion
	double moment_tol;			//Factor applied the norm of the vector of moments of the system for estabelecer convergence criterion
	
	double force_ref;			//Value of reference for evaluation of the stopping criterion of forces
	int force_n_times;			//Number of values used in the calculation of the media of the value of reference
	double moment_ref;			//Value of reference for evaluation of the stopping criterion of moments
	int moment_n_times;			//Number of values used in the calculation of the media of the value of reference


	double force_min;			//Value of force considered small - used when there is only forces nulas - value dimensional
	double moment_min;			//Value of moment considered small - used when there is only moments null - value dimensional
	double constraint_min;		//Value of error of constraint, considered small - always and used when there is joints

	double disp_tol;			//Factor applied the norm of the displacement vector incrementais of the system for estabelecer convergence criterion
	double rot_tol;				//Factor applied the norm of the vector of rotations incrementais of the system for estabelecer convergence criterion
	double lag_tol;				//Factor applied the norm of the vector of Lagrange multipliers of the system for estabelecer convergence criterion

	double disp_min;			//Value of displacement considered small - used when there is only displacements null - value dimensional
	double rot_min;				//Value of rotation considered small - used when there is only rotations null - value dimensional
	double lag_min;				//Value of Lagrange multiplier considered small - used when there is only Lagrange multipliers null - value dimensional

	double divergence_ref;		//Value bastante elevado that indicates divergence of the modelo
	bool diverged;				//Variable booleada that indicates divergence

	double disp_criterion;		//Convergence criterion for norm of displacement - calculated automaticamente
	double rot_criterion;		//Convergence criterion for norm of rotation - calculated automaticamente
	double lag_criterion;		//Convergence criterion for norm of Lagrange multiplier - calculated automaticamente
	double force_criterion;		//Convergence criterion for norm of force - calculated automaticamente
	double moment_criterion;	//Convergence criterion for norm of moment - calculated automaticamente	
	double constraint_criterion;//Convergence criterion for norm of constraint - calculated automaticamente

	bool NaNDetector(double &value);

	void PlotDivergenceReport();//Prints in the tela information for auxiliary the because of the divergence

	//Nodes or SpecialConstraint associated to the maximum error
	int super_node_disp;
	int node_disp;
	int node_rot;
	int sc_lag;
	//Nodes or SpecialConstraint associated to the maximum error
	int super_node_force;
	int node_force;
	int node_moment;
	int sc_constraint;

	int n_conv_evaluated;	//Contador of the number of times that the criterion was estabelecido.

};

