#include "LagrangeSave.h"

#include "Matrix.h"
#include"Database.h"
//Global variables
extern
Database db;

LagrangeSave::LagrangeSave(void)
{
	//Variables saved
	Q_i = new Matrix*[2];					//Tensor rotation in the instant i
	Q_i[0] = new Matrix(3, 3);
	Q_i[1] = new Matrix(3, 3);

	alpha_i = new Matrix*[2];				//Vector alpha in the instant i 
	alpha_i[0] = new Matrix(3);
	alpha_i[1] = new Matrix(3);

	d_alpha_i = new Matrix*[2];				//Vector d_alpha in the instant i 
	d_alpha_i[0] = new Matrix(3);
	d_alpha_i[1] = new Matrix(3);

	theta_i = new Matrix*[2];				//Vector theta in the instant i 
	theta_i[0] = new Matrix(3);
	theta_i[1] = new Matrix(3);

	kappa_i_ref = new Matrix*[2];			//Vector curvature generalizada in the instant i
	kappa_i_ref[0] = new Matrix(3);
	kappa_i_ref[1] = new Matrix(3);

	u_i = new Matrix*[2];					//Vector position of the axis of the bar in the instant i
	u_i[0] = new Matrix(3);
	u_i[1] = new Matrix(3);

	dz_i = new Matrix*[2];					//Vector derivative of the displacement in the instant i
	dz_i[0] = new Matrix(3);
	dz_i[1] = new Matrix(3);

	//Initialization of values in the variables
	(*(Q_i[0]))(0, 0) = 1.0;
	(*(Q_i[0]))(1, 1) = 1.0;
	(*(Q_i[0]))(2, 2) = 1.0;

	(*(Q_i[1]))(0, 0) = 1.0;
	(*(Q_i[1]))(1, 1) = 1.0;
	(*(Q_i[1]))(2, 2) = 1.0;

	(*dz_i[0])(2, 0) = 1.0;		//Always will be this value (system local)
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

//Assignment operator	
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
