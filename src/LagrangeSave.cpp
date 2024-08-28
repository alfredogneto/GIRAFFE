#include "LagrangeSave.h"

#include "Matrix.h"
#include"Database.h"
//Variaveis globais
extern
Database db;

LagrangeSave::LagrangeSave(void)
{
	//Variaveis salvas
	Q_i = new Matrix*[2];					//Tensor rotação no instante i
	Q_i[0] = new Matrix(3, 3);
	Q_i[1] = new Matrix(3, 3);

	alpha_i = new Matrix*[2];				//Vetor alpha no instante i 
	alpha_i[0] = new Matrix(3);
	alpha_i[1] = new Matrix(3);

	d_alpha_i = new Matrix*[2];				//Vetor d_alpha no instante i 
	d_alpha_i[0] = new Matrix(3);
	d_alpha_i[1] = new Matrix(3);

	theta_i = new Matrix*[2];				//Vetor theta no instante i 
	theta_i[0] = new Matrix(3);
	theta_i[1] = new Matrix(3);

	kappa_i_ref = new Matrix*[2];			//Vetor curvatura generalizada no instante i
	kappa_i_ref[0] = new Matrix(3);
	kappa_i_ref[1] = new Matrix(3);

	u_i = new Matrix*[2];					//Vetor posição do eixo da barra no instante i
	u_i[0] = new Matrix(3);
	u_i[1] = new Matrix(3);

	dz_i = new Matrix*[2];					//Vetor derivada do deslocamento no instante i
	dz_i[0] = new Matrix(3);
	dz_i[1] = new Matrix(3);

	//Inicialização de valores nas variaveis
	(*(Q_i[0]))(0, 0) = 1.0;
	(*(Q_i[0]))(1, 1) = 1.0;
	(*(Q_i[0]))(2, 2) = 1.0;

	(*(Q_i[1]))(0, 0) = 1.0;
	(*(Q_i[1]))(1, 1) = 1.0;
	(*(Q_i[1]))(2, 2) = 1.0;

	(*dz_i[0])(2, 0) = 1.0;		//Sempre sera esse valor (sistema local)
	(*dz_i[1])(2, 0) = 1.0;

	(*d_alpha_i[0])(0, 0) = 0;
	(*d_alpha_i[0])(1, 0) = 0;
	(*d_alpha_i[0])(2, 0) = 0;

	(*d_alpha_i[1])(0, 0) = 0;
	(*d_alpha_i[1])(1, 0) = 0;
	(*d_alpha_i[1])(2, 0) = 0;

	(*alpha_i[0])(0, 0) = 0;
	(*alpha_i[0])(1, 0) = 0;
	(*alpha_i[0])(2, 0) = 0;

	(*alpha_i[1])(0, 0) = 0;
	(*alpha_i[1])(1, 0) = 0;
	(*alpha_i[1])(2, 0) = 0;

}
LagrangeSave::~LagrangeSave(void)
{
	delete Q_i[0];
	delete Q_i[1];
	delete[] Q_i;

	delete alpha_i[0];
	delete alpha_i[1];
	delete[] alpha_i;

	delete d_alpha_i[0];
	delete d_alpha_i[1];
	delete[] d_alpha_i;

	delete theta_i[0];
	delete theta_i[1];
	delete[] theta_i;

	delete kappa_i_ref[0];
	delete kappa_i_ref[1];
	delete[] kappa_i_ref;

	delete u_i[0];
	delete u_i[1];
	delete[] u_i;

	delete dz_i[0];
	delete dz_i[1];
	delete[] dz_i;

}

//Operador de Atribuição	
LagrangeSave &LagrangeSave::operator = (LagrangeSave const &save1)
{
	*Q_i[0] = *save1.Q_i[0];
	*Q_i[1] = *save1.Q_i[1];

	*kappa_i_ref[0] = *save1.kappa_i_ref[0];
	*kappa_i_ref[1] = *save1.kappa_i_ref[1];

	*alpha_i[0] = *save1.alpha_i[0];
	*alpha_i[1] = *save1.alpha_i[1];

	*d_alpha_i[0] = *save1.d_alpha_i[0];
	*d_alpha_i[1] = *save1.d_alpha_i[1];

	*theta_i[0] = *save1.theta_i[0];
	*theta_i[1] = *save1.theta_i[1];

	*u_i[0] = *save1.u_i[0];
	*u_i[1] = *save1.u_i[1];

	*dz_i[0] = *save1.dz_i[0];
	*dz_i[1] = *save1.dz_i[1];

	return *this;
}
