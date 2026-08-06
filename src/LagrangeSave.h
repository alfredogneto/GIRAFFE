#pragma once
class Matrix;

class LagrangeSave
{
public:
	LagrangeSave(void);
	~LagrangeSave(void);

	//Variables saved in the gauss points of the elements of viga (for to be used posteriormente without necessidade of interpolation
	Matrix** Q_i;					//Tensor rotation in the instant i
	Matrix** alpha_i;				//Vector rotation in the instant i 
	Matrix** d_alpha_i;				//Derivative temporal of the vector rotation in the instant i
	Matrix** theta_i;				//Vector rotation of Euler in the instant i
	Matrix** kappa_i_ref;			//Vector curvature generalizada in the instant i
	Matrix** u_i;					//Vector displacement of the axis of the bar in the instant i
	Matrix** dz_i;					//Vector derivative of the displacement in the instant i
	
	LagrangeSave &operator = (LagrangeSave const &save1);
};
