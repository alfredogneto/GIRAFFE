#include "Node.h"

#include"Database.h"
//Global variables
extern
Database db;

Node::Node(int e_nGL)
{
	nGL = e_nGL;						//Sets the number of degrees of freedom by node
	ref_coordinates = new double[nGL];	//Coordinates of the node in the configuration of reference
	copy_coordinates = new double[nGL];	//Coordinates of the node in the configuration of copies (last converged)
	copy_rot_euler = new double[3];
	displacements = new double[nGL];	//Displacements of the node in relação the last copies of coordinates
	vel = new double[nGL];				//Velocities
	copy_vel = new double[nGL];		
	accel = new double[nGL];			//Acelerações
	copy_accel = new double[nGL];

	constraints = new int[nGL];
	GLs = new int[nGL];
	active_GL = new int[nGL];

	//Zeroing values initial
	for (int i = 0; i < nGL; i++)
	{
		GLs[i] = 0;
		displacements[i] = 0;
		ref_coordinates[i] = 0;
		copy_coordinates[i] = 0;
		vel[i] = 0;
		accel[i] = 0;
		copy_vel[i] = 0;
		copy_accel[i] = 0;
		active_GL[i] = 0;
		constraints[i] = 0;
	}
	copy_rot_euler[0] = 0.0;
	copy_rot_euler[1] = 0.0;
	copy_rot_euler[2] = 0.0;

	n_GL_free = 0;
	n_GL_fixed = 0;

	theta = 0;
	rot_euler = new Matrix(3);
	rot_rodrigues = new Matrix(3);
	load = new Matrix(db.number_GLs_node);
	moment = new Matrix(3);
	force = new Matrix(3);

	A = new Matrix(3,3);
	g = 0;
	I = new Matrix(3,3);
	(*I)(0, 0) = 1.0;
	(*I)(1, 1) = 1.0;
	(*I)(2, 2) = 1.0;
	Xi = new Matrix (3,3);
	Xi_T_inv = new Matrix(3,3);
	Q = new Matrix(3, 3);
	Q0 = new Matrix(3, 3);
	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;
	//Initialization of the matrix Q
	*Q = *I3;
	//Initialization of the matrix Q0
	*Q0 = *I3;

	alpha_escalar = 0;
	theta_escalar = 0;
	alpha_v = new Matrix(3, 1);
	theta_v = new Matrix(3, 1);

	alpha_1 = Matrix(3);
	alpha_2 = Matrix(3);
	alpha_3 = Matrix(3);

	theta3 = Matrix(3);
	
	//Explicit
	u = new Matrix(3, 1);
	du = new Matrix(3, 1);
	ddu = new Matrix(3, 1);
	alpha = new Matrix(3, 1);
	omega = new Matrix(3, 1);
	domega = new Matrix(3, 1);
	
	copy_du = new Matrix(3, 1);
	copy_omega = new Matrix(3, 1);

	flag_material_description = false;
	flag_pseudo_moment = true;
}

void Node::ZeroIC()
{
	//Zeroing values initial
	for (int i = 0; i < nGL; i++)
	{
		vel[i] = 0;
		accel[i] = 0;
		copy_vel[i] = 0;
		copy_accel[i] = 0;
	}
}


Node::~Node()
{
	delete[] ref_coordinates;
	delete[] copy_coordinates;
	delete[] copy_rot_euler;
	delete[] displacements;
	delete[] vel;
	delete[] copy_vel;
	delete[] accel;
	delete[] copy_accel;

	delete[] constraints;
	delete[] GLs;
	delete[] active_GL;

	delete rot_euler;
	delete rot_rodrigues;
	delete load;
	delete moment;
	delete force;
	delete A;
	delete I;
	delete Xi;
	delete Xi_T_inv;
	delete Q;
	delete I3;
	delete Q0;

	delete alpha_v;
	delete theta_v;

	delete u;
	delete du;
	delete ddu;
	delete alpha;
	delete omega;
	delete domega;
	
	delete copy_du;
	delete copy_omega;
}

bool Node::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Checks the keyword "Node"
	if (!strcmp(s, "Node"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;

	//Reading of the coordinates
	fscanf(f, "%s", s);
	ref_coordinates[0] = atof(s);
	copy_coordinates[0] = atof(s);

	fscanf(f, "%s", s);
	ref_coordinates[1] = atof(s);
	copy_coordinates[1] = atof(s);

	fscanf(f, "%s", s);
	ref_coordinates[2] = atof(s);
	copy_coordinates[2] = atof(s);

	return true; 
}

void Node::Write(FILE *f)
{
	fprintf(f, "Node\t%d\t%.12e\t%.12e\t%.12e\n", number,
		ref_coordinates[0],
		ref_coordinates[1],
		ref_coordinates[2]);
}

void Node::WriteUpdated(FILE *f)
{
	fprintf(f, "Node\t%d\t%.12e\t%.12e\t%.12e\n", number,
		copy_coordinates[0],
		copy_coordinates[1],
		copy_coordinates[2]);
}

void Node::WriteInitialConditions(FILE *f)
{
	fprintf(f, "Node\t%d\tDU\t%.12e\t%.12e\t%.12e\tOMEGA\t%.12e\t%.12e\t%.12e\t", number,
		copy_vel[0],
		copy_vel[1],
		copy_vel[2],
		copy_vel[3],
		copy_vel[4],
		copy_vel[5]);
}

void Node::WriteResults(FILE *f)
{
	fprintf(f, "Node\t%d\tX\t%.6e\tY\t%.6e\tZ\t%.6e\tRX\t%.6e\tRY\t%.6e\tRZ\t%.6e\n", number,
	copy_coordinates[0],
	copy_coordinates[1],
	copy_coordinates[2],
	copy_rot_euler[0],
	copy_rot_euler[1],
	copy_rot_euler[2]);

	//Acrescentar here other degrees of freedom possiveis
}

void Node::WriteVTK(FILE *f)
{
	fprintf(f, "%lf %lf %lf\n",
		copy_coordinates[0],
		copy_coordinates[1],
		copy_coordinates[2]);
}
void Node::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Calculation of the angle of rotation
	(*rot_rodrigues)(0, 0) = displacements[3];
	(*rot_rodrigues)(1, 0) = displacements[4];
	(*rot_rodrigues)(2, 0) = displacements[5];
	//Internal forces
	for (int j = 0; j < db.number_GLs_node; j++)
	{
		if (GLs[j] < 0 && active_GL[j] == 1)	//If the degree of freedom for fixed and active
			(*load)(j, 0) = db.global_P_B(-GLs[j] - 1, 0);
		if (GLs[j] > 0 && active_GL[j] == 1)	//If the degree of freedom for free and active
			(*load)(j, 0) = db.global_I_A(+GLs[j] - 1, 0);
	}
	
	Matrix local_vel(3);
	Matrix local_ang_vel(3);
	Matrix local_accel(3);
	Matrix local_ang_accel(3);

	for (int j = 0; j < 3; j++)
	{
		(*force)(j, 0) = (*load)(j, 0);
		(*moment)(j, 0) = (*load)(j + 3, 0);
	}
		
	if (flag_pseudo_moment)
	{
		//Conversion of the pseudo-moment for moment (parameters of rotation of Rodrigues)
		//Calculating the operator Xi
		*A = skew(*rot_rodrigues);			//Matrix A
		g = 4.0 / (4.0 + norm(*rot_rodrigues)*norm(*rot_rodrigues));		//function g(alpha) - in some ref. tb. called of h(alpha)
		*Xi = g * (*I + 0.5*(*A));
		*Xi_T_inv = invert3x3(transp(*Xi));
		*moment = (*Xi_T_inv)*(*moment);
		for (int j = 0; j < 3; j++)
			(*load)(j + 3, 0) = (*moment)(j, 0);
	}
	
	Matrix global_force;
	Matrix global_moment;
	if (flag_material_description)
	{
		global_force = (*Q)*(*Q0)*(*force);
		global_moment = (*Q)*(*Q0)*(*moment);
	}
	else
	{
		global_force = *force;
		global_moment = *moment;
	}
		

	//Header
	if (first_record == true)
		fprintf(f, "TIME\tX\tY\tZ\tROTX\tROTY\tROTZ\tFX\tFY\tFZ\tMX\tMY\tMZ\tFXloc\tFYloc\tFZloc\tMXloc\tMYloc\tMZloc\tdX\tdY\tdZ\tdROTX\tdROTY\tdROTZ\tddX\tddY\tddZ\tddROTX\tddROTY\tddROTZ\n");
	//Information the to be saved
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
		time,
		//Displacements - global, rotations can be local or global
		copy_coordinates[0],
		copy_coordinates[1],
		copy_coordinates[2],
		copy_rot_euler[0],
		copy_rot_euler[1],
		copy_rot_euler[2],
		//Forces - global
		global_force(0, 0),
		global_force(1, 0),
		global_force(2, 0),
		//Moments - global
		global_moment(0, 0),
		global_moment(1, 0),
		global_moment(2, 0),
		//Forces and moments - local
		(*load)(0, 0),
		(*load)(1, 0),
		(*load)(2, 0),
		(*load)(3, 0),
		(*load)(4, 0),
		(*load)(5, 0),
		copy_vel[0],
		copy_vel[1],
		copy_vel[2],
		copy_vel[3],
		copy_vel[4],
		copy_vel[5],

		copy_accel[0],
		copy_accel[1],
		copy_accel[2],
		copy_accel[3],
		copy_accel[4],
		copy_accel[5]
		);
}

void Node::SaveConfiguration()
{
	//Update of displacements (description lagrangiana updated - incremental)
	copy_coordinates[0] += displacements[0];
	copy_coordinates[1] += displacements[1];
	copy_coordinates[2] += displacements[2];

	//Update of rotations (description lagrangiana updated - incremental)
	alpha_1(0, 0) = copy_coordinates[3];
	alpha_1(1, 0) = copy_coordinates[4];
	alpha_1(2, 0) = copy_coordinates[5];
	alpha_2(0, 0) = displacements[3];
	alpha_2(1, 0) = displacements[4];
	alpha_2(2, 0) = displacements[5];

	//Fórmula of Rodrigues
	if (flag_material_description == false)
		alpha_3 = 4.0 / (4.0 - dot(alpha_2, alpha_1))*(alpha_2 + alpha_1 + 0.5*cross(alpha_2, alpha_1));	//Vector rotation of Rodrigues - update of description espacial
	else
		alpha_3 = 4.0 / (4.0 - dot(alpha_1, alpha_2))*(alpha_1 + alpha_2 + 0.5*cross(alpha_1, alpha_2));	//Vector rotation of Rodrigues - update of description material
			
	copy_coordinates[3] = alpha_3(0, 0);
	copy_coordinates[4] = alpha_3(1, 0);
	copy_coordinates[5] = alpha_3(2, 0);
	

	////////////////////Saving rotations of Euler////////////////////
	theta_escalar3 = 2.0*atan(norm(alpha_3) / 2.0);
	if (norm(alpha_3) != 0.0)
		theta3 = (theta_escalar3 / norm(alpha_3))*(alpha_3);
	else
		zeros(&theta3);

	copy_rot_euler[0] = theta3(0, 0);
	copy_rot_euler[1] = theta3(1, 0);
	copy_rot_euler[2] = theta3(2, 0);

	////////////////////////////////////////////////////////////////////////////////////////
	alpha_escalar = norm(alpha_3);
	*A = skew(alpha_3);
	g = 4.0 / (4.0 + alpha_escalar * alpha_escalar);
	*Q = *I + g * (*A + 0.5*(*A)*(*A));
	////////////////////////////////////////////////////////////////////////////////////////
	
	//Update of other degrees of freedom (description lagrangiana total)
	for (int j = 6; j < db.number_GLs_node; j++)
	{
		copy_coordinates[j] = displacements[j];
	}

	//Update of the velocities and acelerações
	for (int j = 0; j < db.number_GLs_node; j++)
	{
		copy_accel[j] = accel[j];
		copy_vel[j] = vel[j];
	}
	//Explicit
	*copy_du = (*du);
	*copy_omega = *omega;
}

Matrix Node::InvXi(const Matrix& alpha)
{
	double v[50];
	Matrix invXi(3, 3);

	if (flag_material_description)
	{
		v[1] = alpha.m_matrix[0];
		v[2] = alpha.m_matrix[1];
		v[3] = alpha.m_matrix[2];
		v[13] = 4e0 / (4e0 + (v[1] * v[1]) + (v[2] * v[2]) + (v[3] * v[3]));
		v[35] = 0.5e0*v[13];
		v[18] = (v[13] * v[13]);
		v[14] = v[3] * v[35];
		v[22] = -(v[13] * v[14]);
		v[15] = -(v[2] * v[35]);
		v[30] = -(v[14] * v[15]);
		v[27] = -(v[13] * v[15]);
		v[16] = 1e0*v[1] * v[35];
		v[31] = -(v[13] * v[16]);
		v[28] = v[14] * v[16];
		v[23] = -(v[15] * v[16]);
		invXi(0, 0) = ((v[16] * v[16]) + v[18]) / v[18];
		invXi(0, 1) = (v[22] + v[23]) / v[18];
		invXi(0, 2) = (v[27] + v[28]) / v[18];
		invXi(1, 0) = (-v[22] + v[23]) / v[18];
		invXi(1, 1) = ((v[15] * v[15]) + v[18]) / v[18];
		invXi(1, 2) = (v[30] + v[31]) / v[18];
		invXi(2, 0) = (-v[27] + v[28]) / v[18];
		invXi(2, 1) = (v[30] - v[31]) / v[18];
		invXi(2, 2) = ((v[14] * v[14]) + v[18]) / v[18];
	}
	else
	{
		v[1] = alpha.m_matrix[0];
		v[2] = alpha.m_matrix[1];
		v[3] = alpha.m_matrix[2];
		v[13] = 4e0 / (4e0 + (v[1] * v[1]) + (v[2] * v[2]) + (v[3] * v[3]));
		v[26] = 0.5e0*v[13];
		v[18] = (v[13] * v[13]);
		v[14] = -(v[26] * v[3]);
		v[15] = v[2] * v[26];
		v[16] = -1e0*v[1] * v[26];
		v[19] = v[13] * v[14];
		v[20] = -(v[15] * v[16]);
		v[21] = v[13] * v[15];
		v[22] = v[14] * v[16];
		v[23] = -(v[14] * v[15]);
		v[24] = v[13] * v[16];
		invXi(0, 0) = ((v[16] * v[16]) + v[18]) / v[18];
		invXi(0, 1) = (-v[19] + v[20]) / v[18];
		invXi(0, 2) = (-v[21] + v[22]) / v[18];
		invXi(1, 0) = (v[19] + v[20]) / v[18];
		invXi(1, 1) = ((v[15] * v[15]) + v[18]) / v[18];
		invXi(1, 2) = (v[23] - v[24]) / v[18];
		invXi(2, 0) = (v[21] + v[22]) / v[18];
		invXi(2, 1) = (v[23] + v[24]) / v[18];
		invXi(2, 2) = ((v[14] * v[14]) + v[18]) / v[18];
	}
	
	return invXi;
}