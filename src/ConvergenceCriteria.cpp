#pragma once
#include "ConvergenceCriteria.h"

#include "Node.h"
#include "SuperNode.h"
#include "SpecialConstraint.h"

#include"Database.h"
//Variaveis globais
extern
Database db;



ConvergenceCriteria::ConvergenceCriteria()
{
	force_tol = 1e-4;			//Fator aplicado a norma do vetor de forças do sistema para estabelecer criterio de convergência
	moment_tol = 1e-4;			//Fator aplicado a norma do vetor de momentos do sistema para estabelecer criterio de convergência

	force_min = 1e-5;			//Valor de força considerado pequeno - usado quando ha somente forças nulas - valor dimensional
	moment_min = 1e-5;			//Valor de momento considerado pequeno - usado quando ha somente momentos nulos - valor dimensional
	constraint_min = 1e-7;		//Valor de erro de constraint, considerado pequeno - sempre e utilizado quando ha joints

	disp_tol = 1e-4;			//Fator aplicado a norma do vetor de deslocamentos incrementais do sistema para estabelecer criterio de convergência
	rot_tol = 1e-4;				//Fator aplicado a norma do vetor de rotações incrementais do sistema para estabelecer criterio de convergência
	lag_tol = 1e-4;				//Fator aplicado a norma do vetor de multiplicadores de lagrange do sistema para estabelecer criterio de convergência

	disp_min = 1e-6;			//Valor de deslocamento considerado pequeno - usado quando ha somente deslocamentos nulos - valor dimensional
	rot_min = 1e-6;				//Valor de rotação considerado pequeno - usado quando ha somente rotações nulos - valor dimensional
	lag_min = 1e-6;				//Valor de multiplicador de lagrange considerado pequeno - usado quando ha somente multiplicadores de lagrange nulos - valor dimensional

	divergence_ref = 1e15;		//Valor bastante elevado que indica divergência do modelo

	disp_criterion = 0.0;		//Criterio de convergência para norma de deslocamento - calculado automaticamente
	rot_criterion = 0.0;		//Criterio de convergência para norma de rotação - calculado automaticamente
	lag_criterion = 0.0;		//Criterio de convergência para norma de multiplicador de lagrange - calculado automaticamente
	force_criterion = 0.0;		//Criterio de convergência para norma de força - calculado automaticamente
	moment_criterion = 0.0;		//Criterio de convergência para norma de momento - calculado automaticamente	
	constraint_criterion = 0.0;	//Criterio de convergência para norma de constraint - calculado automaticamente

	diverged = false;

	node_disp = 0;
	super_node_disp = 0;
	node_rot = 0;
	sc_lag = 0;
	node_force = 0;
	super_node_force = 0;
	node_moment = 0;
	sc_constraint = 0;

	force_ref = 0.0;			//Valor de referência para avaliação do criterio de parada de forças
	force_n_times = 0;			//Numero de valores utilizado no calculo da media do valor de referência
	moment_ref = 0.0;			//Valor de referência para avaliação do criterio de parada de momentos
	moment_n_times = 0;			//Numero de valores utilizado no calculo da media do valor de referência

	n_conv_evaluated = 0;
}

ConvergenceCriteria::~ConvergenceCriteria()
{
}

bool ConvergenceCriteria::Read(FILE *f)
{
	char s[100];
	fscanf(f, "%s", s);
	if (!strcmp(s, "ForceTolerance"))
	{
		fscanf(f, "%s", s);
		force_tol = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "MomentTolerance"))
	{
		fscanf(f, "%s", s);
		moment_tol = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "ForceMinimumReference"))
	{
		fscanf(f, "%s", s);
		force_min = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "MomentMinimumReference"))
	{
		fscanf(f, "%s", s);
		moment_min = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "ConstraintMinimumReference"))
	{
		fscanf(f, "%s", s);
		constraint_min = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "DisplacementTolerance"))
	{
		fscanf(f, "%s", s);
		disp_tol = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "RotationTolerance"))
	{
		fscanf(f, "%s", s);
		rot_tol = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "LagrangeTolerance"))
	{
		fscanf(f, "%s", s);
		lag_tol = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "DisplacementMinimumReference"))
	{
		fscanf(f, "%s", s);
		disp_min = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "RotationMinimumReference"))
	{
		fscanf(f, "%s", s);
		rot_min = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "LagrangeMinimumReference"))
	{
		fscanf(f, "%s", s);
		lag_min = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "DivergenceReference"))
	{
		fscanf(f, "%s", s);
		divergence_ref = atof(s);
	}
	else
		return false;

	return true;
}

void ConvergenceCriteria::Write(FILE *f)
{
	//Escrita
	fprintf(f, "ForceTolerance\t%.6e\n",force_tol);
	fprintf(f, "MomentTolerance\t%.6e\n", moment_tol);

	fprintf(f, "ForceMinimumReference\t%.6e\n", force_min);
	fprintf(f, "MomentMinimumReference\t%.6e\n", moment_min);
	fprintf(f, "ConstraintMinimumReference\t%.6e\n", constraint_min);

	fprintf(f, "DisplacementTolerance\t%.6e\n", disp_tol);
	fprintf(f, "RotationTolerance\t%.6e\n", rot_tol);
	fprintf(f, "LagrangeTolerance\t%.6e\n", lag_tol);

	fprintf(f, "DisplacementMinimumReference\t%.6e\n", disp_min);
	fprintf(f, "RotationMinimumReference\t%.6e\n", rot_min);
	fprintf(f, "LagrangeMinimumReference\t%.6e\n", lag_min);

	fprintf(f, "DivergenceReference\t%.6e\n", divergence_ref);
}
void ConvergenceCriteria::EstablishResidualCriteria()
{
	//Maximos incrementos
	double max_force = 0;
	double max_moment = 0;

	/*double sum_forces = 0.0;
	double sum_moments = 0.0;
	int n_forces = 0;
	int n_moments = 0;*/


	int GL = 0;
	double value = 0;//variavel temporaria para salvar valores de GL
	//Varredura dos nós
	for (int node = 0; node < db.number_nodes; node++)
	{
		for (int i = 0; i < db.number_GLs_node; i++)
		{
			GL = db.nodes[node]->GLs[i];
			if (GL > 0)
			{
				value = abs(db.global_P_A(GL - 1, 0));
				//displacement DOF
				if (i == 0 || i == 1 || i == 2)
				if (value > max_force)
					max_force = value;
				//rotation DOF
				if (i == 3 || i == 4 || i == 5)
				if (value > max_moment)
					max_moment = value;	
				////New criterion
				//value = abs(db.global_ABS_P_A(GL - 1, 0));
				////displacement DOF
				//if (i == 0 || i == 1 || i == 2)
				//{
				//	sum_forces += value;
				//	n_forces += db.global_COUNT_ABS_P_A(GL - 1, 0);
				//}
				////rotation DOF
				//if (i == 3 || i == 4 || i == 5)
				//{
				//	sum_moments += value;
				//	n_moments += db.global_COUNT_ABS_P_A(GL - 1, 0);
				//}

			}
		}
	}
	//Varredura dos super nodes
	for (int s_node = 0; s_node < db.number_super_nodes; s_node++)
	{
		for (int i = 0; i < db.super_nodes[s_node]->n_displacement_DOFs; i++)
		{
			GL = db.super_nodes[s_node]->DOFs[i];
			if (GL > 0)
			{
				value = abs(db.global_P_A(GL - 1, 0));
				//displacement DOF
				if (value > max_force)
					max_force = value;
				////New criterion
				//value = abs(db.global_ABS_P_A(GL - 1, 0));
				//sum_forces += value;
				//n_forces += db.global_COUNT_ABS_P_A(GL - 1, 0);
			}
		}
	}

	double prev_force_criterion = force_criterion;
	double prev_moment_criterion = moment_criterion;

	//Estabelecendo os criterios
	if (max_force*force_tol > force_min*force_tol)
		force_criterion = max_force*force_tol;
	else
		force_criterion = force_min*force_tol;
	if (max_moment*moment_tol > moment_min*moment_tol)
		moment_criterion = max_moment*moment_tol;
	else
		moment_criterion = moment_min*moment_tol;
	constraint_criterion = constraint_min;

	//Updating criteria: historical data
	force_criterion = (n_conv_evaluated * prev_force_criterion + force_criterion) / (n_conv_evaluated + 1);
	moment_criterion = (n_conv_evaluated * prev_moment_criterion + moment_criterion) / (n_conv_evaluated + 1);

	//double max_force_new = sum_forces / n_forces;
	//double max_moment_new = sum_moments / n_moments;

	////Estabelecendo os criterios - new
	//if (max_force_new*force_tol > force_min*force_tol)
	//	force_criterion = max_force_new * force_tol;
	//else
	//	force_criterion = force_min * force_tol;
	//if (max_moment_new*moment_tol > moment_min*moment_tol)
	//	moment_criterion = max_moment_new * moment_tol;
	//else
	//	moment_criterion = moment_min * moment_tol;
	
	n_conv_evaluated++;
}
bool ConvergenceCriteria::CheckGLConvergence()
{
	bool converged = true;
	int GL = 0;
	//Maximos incrementos
	double max_disp_inc = 0;
	double max_rot_inc = 0;
	double max_lag_inc = 0;
	//Valores dos GLs
	double max_disp_value = 0;
	double max_rot_value = 0;
	double max_lag_value = 0;


	bool plot_force_disp = false;
	bool plot_moment_rot = false;
	
	double inc = 0;//variavel temporaria para salvar incrementos de GL
	double value = 0;//variavel temporaria para salvar valores de GL
	//Varredura dos nós
	for (int node = 0; node < db.number_nodes; node++)
	{
		for (int i = 0; i < db.number_GLs_node; i++)
		{
			GL = db.nodes[node]->GLs[i];
			if (GL > 0)
			{
				inc = abs(db.global_P_A(GL-1, 0));
				if (NaNDetector(inc))
					diverged = true;
				value = abs(db.nodes[node]->displacements[i]);
				//displacement
				if (i == 0 || i == 1 || i == 2)
				{ 
					plot_force_disp = true;
					if (inc > max_disp_inc)
					{
						max_disp_inc = inc;
						node_disp = node + 1;
						super_node_disp = 0;
					}
					if (value > max_disp_value)
						max_disp_value = value;
				}
				//rotation
				if (i == 3 || i == 4 || i == 5)
				{
					plot_moment_rot = true;
					if (inc > max_rot_inc)
					{
						max_rot_inc = inc;
						node_rot = node + 1;
					}
					if (value > max_rot_value)
						max_rot_value = value;
				}
			}
		}
	}

	//Varredura dos super nodes
	for (int s_node = 0; s_node < db.number_super_nodes; s_node++)
	{
		for (int i = 0; i < db.super_nodes[s_node]->n_displacement_DOFs; i++)
		{
			GL = db.super_nodes[s_node]->DOFs[i];
			if (GL > 0)
			{
				plot_force_disp = true;
				inc = abs(db.global_P_A(GL - 1, 0));
				if (NaNDetector(inc))
					diverged = true;
				value = abs(db.super_nodes[s_node]->displacements[i]);
				//displacement
				if (inc > max_disp_inc)
				{
					max_disp_inc = inc;
					super_node_disp = s_node + 1;
					node_disp = 0;
				}
				if (value > max_disp_value)
					max_disp_value = value;
			}
		}
	}

	//varredura dos special constraints
	if (db.special_constraints_exist == true)
	{
		for (int sc = 0; sc < db.number_special_constraints; sc++)
		{
			for (int i = 0; i < db.special_constraints[sc]->n_GL; i++)
			{
				GL = db.special_constraints[sc]->GLs[i];
				if (GL > 0)
				{
					inc = abs(db.global_P_A(GL - 1, 0));
					if (NaNDetector(inc))
						diverged = true;
					value = abs(db.special_constraints[sc]->lambda[i]);
					//special constraint
					if (inc > max_lag_inc)
					{
						max_lag_inc = inc;
						sc_lag = sc + 1;
					}
					if (value > max_lag_value)
						max_lag_value = value;
				}
			}
		}
	}
	
	
	//Estabelecendo os criterios de convergência
	//Deslocamentos
	if (max_disp_value*disp_tol > disp_min*disp_tol)
		disp_criterion = disp_tol*max_disp_value;
	else
		disp_criterion = disp_min*disp_tol;
	//Rotações
	if (max_rot_value*rot_tol > rot_min*rot_tol)
		rot_criterion = rot_tol*max_rot_value;
	else
		rot_criterion = rot_min*rot_tol;
	//Lagrange
	if (max_lag_value*lag_tol > lag_min*lag_tol)
		lag_criterion = lag_tol*max_lag_value;
	else
		lag_criterion = lag_min*lag_tol;

	//Testes de divergência
	if (disp_criterion > divergence_ref || rot_criterion > divergence_ref || lag_criterion > divergence_ref)
	{
		diverged = true;
		converged = false;
		return converged;
	}
		
	if (plot_force_disp == true)
	{
		//Testes de convergência
		db.myprintf("\tDisplacement: Residual %.3e Criterion %.3e  ", max_disp_inc, disp_criterion);
		if (max_disp_inc > disp_criterion)
			converged = false;
		else
			db.myprintf("CONVERGED");
		db.myprintf("\n");
	}
	if (plot_moment_rot == true)
	{
		db.myprintf("\tRotation    : Residual %.3e Criterion %.3e  ", max_rot_inc, rot_criterion);
		if (max_rot_inc > rot_criterion)
			converged = false;
		else
			db.myprintf("CONVERGED");
		db.myprintf("\n");
	}
	
	if (db.special_constraints_exist == true)
	{
		db.myprintf("\tLagrange    : Residual %.3e Criterion %.3e  ", max_lag_inc, lag_criterion);
		if (max_lag_inc > lag_criterion)
			converged = false;
		else
			db.myprintf("CONVERGED");
		db.myprintf("\n");
	}
	return converged;
}

bool ConvergenceCriteria::CheckResidualConvergence()
{
	bool converged = true;
	int GL = 0;
	//Maximos incrementos
	double max_force = 0;
	double max_moment = 0;
	double max_constraint = 0;


	bool plot_force_disp = false;
	bool plot_moment_rot = false;
	
	double value = 0;//variavel temporaria para salvar valores de GL
	double error = 0;
	//Varredura dos nós
	for (int node = 0; node < db.number_nodes; node++)
	{
		for (int i = 0; i < db.number_GLs_node; i++)
		{
			GL = db.nodes[node]->GLs[i];
			if (GL > 0)
			{
				value = abs(db.global_P_A(GL - 1, 0));
				if (NaNDetector(value))
					diverged = true;
				//displacement
				if (i == 0 || i == 1 || i == 2)
				{
					plot_force_disp = true;
					if (value > max_force)
					{
						max_force = value;
						node_force = node + 1;
						super_node_force = 0;
					}
				}
				if (i == 3 || i == 4 || i == 5)//rot
				{
					plot_moment_rot = true;
					if (value > max_moment)
					{
						max_moment = value;
						node_moment = node + 1;
					}
				}
			}
		}
	}
	//Varredura dos super nodes
	for (int s_node = 0; s_node < db.number_super_nodes; s_node++)
	{
		for (int i = 0; i < db.super_nodes[s_node]->n_displacement_DOFs; i++)
		{
			GL = db.super_nodes[s_node]->DOFs[i];
			if (GL > 0)
			{
				plot_force_disp = true;
				value = abs(db.global_P_A(GL - 1, 0));
				if (NaNDetector(value))
					diverged = true;
				//displacement
				if (value > max_force)
				{
					max_force = value;
					super_node_force = s_node + 1;
					node_force = 0;
				}
			}
		}
	}
	//varredura dos special constraints
	if (db.special_constraints_exist == true)
	{
		for (int sc = 0; sc < db.number_special_constraints; sc++)
		{
			for (int i = 0; i < db.special_constraints[sc]->n_GL; i++)
			{
				GL = db.special_constraints[sc]->GLs[i];
				if (GL > 0)
				{
					value = abs(db.global_P_A(GL - 1, 0));
					if (NaNDetector(value))
						diverged = true;
					//special constraint
					if (value > max_constraint)
					{
						max_constraint = value;
						sc_constraint = sc + 1;
					}
				}
			}
		}
	}

	//Testes de divergência
	if (max_force > divergence_ref || max_moment > divergence_ref || max_constraint > divergence_ref)
	{
		diverged = true;
		converged = false;
		return converged;
	}
	
	//Testes de convergência

	if (plot_force_disp == true)
	{
		//Força
		db.myprintf("\tForce       : Residual %.3e Criterion %.3e  ", max_force, force_criterion);
		if (max_force > force_criterion)
			converged = false;
		else
			db.myprintf("CONVERGED");
		db.myprintf("\n");
	}
	if (plot_moment_rot == true)
	{

		//Momento
		db.myprintf("\tMoment      : Residual %.3e Criterion %.3e  ", max_moment, moment_criterion);
		if (max_moment > moment_criterion)
			converged = false;
		else
			db.myprintf("CONVERGED");
		db.myprintf("\n");
	}
	
	if (db.special_constraints_exist == true)
	{
		//Constraint
		db.myprintf("\tConstraint  : Residual %.3e Criterion %.3e  ", max_constraint, constraint_criterion);
		if (max_constraint > constraint_criterion)
			converged = false;
		else
			db.myprintf("CONVERGED");
		db.myprintf("\n");
	}
	return converged;
}

bool ConvergenceCriteria::NaNDetector(double &value)
{
	//Detecção de NaN
	if (value == value)
	{
		//Detecção de valor muito elevado (erro)
		if (value >= 1e300 || value <= -1e300)
			return true;//e NaN ou numero muito elevado
		else
			return false;//Não e NaN
	}
	else
	{
		//db.myprintf("NaN detected!\n");
		return true;//e NaN
	}
		
}
//Imprime na tela informações para auxiliar o porquê da divergência
void ConvergenceCriteria::PlotDivergenceReport()
{
	db.myprintf("\nThe model has diverged! Applying Bisection on loading.\n");
	if (super_node_force == 0)
		db.myprintf("Max Force error occurred at node %d.\n", node_force);
	if (node_force == 0)
		db.myprintf("Max Force error occurred at super node %d.\n", super_node_force);
	db.myprintf("Max Moment error occurred at node %d.\n", node_moment);
	if (db.special_constraints_exist == true)
		db.myprintf("Max Constraint error occurred at special constraint %d.\n", sc_constraint);
	if (super_node_disp == 0)
		db.myprintf("Max Displacement error occurred at node %d.\n", node_disp);
	if (node_disp == 0)
		db.myprintf("Max Displacement error occurred at super node %d.\n", super_node_disp);
	db.myprintf("Max Rotation error occurred at node %d.\n", node_rot);
	if (db.special_constraints_exist == true)
		db.myprintf("Max Lagrange multiplier error occurred at special constraint %d.\n", sc_lag);
}
