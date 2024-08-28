#include "SecHelicalFiber.h"

#include "SolidSection.h"
#include "Hooke.h"
#include "Database.h"
//Variáveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795
SecHelicalFiber::SecHelicalFiber()
{
	sec_details = new SolidSection();

	AC = Matrix(2);
	aero_length = 0;

	R = 0;			//helix radius
	r = 0;			//cross-section radius
	alpha = 0;		//helix angle
	nt = 0;			//number of turns
	nh = 0;			//number of helices
}

SecHelicalFiber::~SecHelicalFiber()
{
	delete sec_details;
}

bool SecHelicalFiber::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "Rh"))
	{
		fscanf(f, "%s", s);
		R = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Rs"))
	{
		fscanf(f, "%s", s);
		r = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Alpha"))
	{
		fscanf(f, "%s", s);
		alpha = atof(s)*PI/180;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nt"))
	{
		fscanf(f, "%s", s);
		nt = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nh"))
	{
		fscanf(f, "%s", s);
		nh = atoi(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "AD"))
	{
		fscanf(f, "%s", s);
		aerodynamicdataID = atoi(s);

		fscanf(f, "%s", s);
		if (!strcmp(s, "AC"))
		{
			fscanf(f, "%s", s);
			AC(0, 0) = atof(s);
			fscanf(f, "%s", s);
			AC(1, 0) = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "AeroLength"))
		{
			fscanf(f, "%s", s);
			aero_length = atof(s);
		}
		else
			return false;
	}
	else
		fsetpos(f, &pos);	//volta à posição anterior

	return true;
}

void SecHelicalFiber::Write(FILE *f)
{
	fprintf(f, "HelicalFiber\t%d\tRh\t%.6e\tRs\t%.6e\tAlpha\t%.6e\tNt\t%d\tNh\t%d\tAD\t%d\n", number, R, r, alpha, nt, nh, aerodynamicdataID);
}

void SecHelicalFiber::PreCalc()
{
	//Section details (for drawing purposes)
	SolidSection* ptr_sd = static_cast<SolidSection*>(sec_details);
	int n_circ = 24;
	ptr_sd->Alloc(n_circ);
	double theta = 0;
	for (int index = 0; index < n_circ; index++)
	{
		theta = (index * 2 * PI) / n_circ;
		ptr_sd->points[index][0] = R*cos(theta);
		ptr_sd->points[index][1] = R*sin(theta);
	}
}

//Fills contents of pointers to D, Mr and Jr (needed for Beam_1)
void SecHelicalFiber::ComputeStiffnessMass(Matrix& D, Matrix& Mr, Matrix& Jr, int material_ID)
{
	Hooke* hooke = static_cast<Hooke*>(db.materials[material_ID - 1]);
	double E = hooke->E;
	double nu = hooke->nu;
	double G = E/(2*(1+nu));
	double rho = hooke->rho;

	// Function - Michele
	// Reference: Marino and Vairo, 2014 :

	// Helix pitch :
	double theta = alpha;
	double p = 2 * PI*R*(1.0/(tan(theta)));
	double c = p / (2 * PI);
	double a = sqrt(R*R + c*c);

	double A = PI*r*r;
	double I1 = PI*r*r*r*r / 4;
	double I2 = I1;
	double J_th = G*PI*r*r*r*r / 2;

	double C3333 = E;
	double S1313 = 1.0/G;
	double S2323 = 1.0 / G;
	double chi11 = 6 * (1 + nu) / (7 + 6 * nu);
	double chi22 = chi11;

	// Compliance matrix(single helix)
	Matrix Ch(6, 6);

	Ch(0, 0) = (PI*nt / a)*(R * R / (C3333*A) + a * a * c * c * (8 * PI * PI * nt * nt - 3) / (6 * C3333*I1)
		+ (c*c*c*c * (8 * PI * PI * nt * nt + 3) + 6 * R * R * (R * R - c * c)) / (6 * C3333*I2)
		+ (a * a * S1313*chi11 + c * c * S2323*chi22) / A
		+ R * R * c * c * (8 * PI * PI * nt * nt + 15) / (6 * J_th));

	Ch(0, 1) = (PI * PI * nt * nt * c * c / a)*(-a * a / (C3333*I1) + c * c / (C3333*I2) + R * R / J_th);
	Ch(1, 0) = Ch(0, 1);

	Ch(0, 2) = -Ch(0, 1)*R / c - 2 * PI * PI * nt * nt * R*c*a / (C3333*I1);
	Ch(2, 0) = Ch(0, 2);

	Ch(0, 3) = (PI*nt*c / (2 * a))*(a * a / (C3333*I1) + (2 * R * R - c * c) / (C3333*I2) - 3 * R * R / J_th);
	Ch(3, 0) = Ch(0, 3);

	Ch(0, 4) = -Ch(0, 2) / R;
	Ch(4, 0) = Ch(0, 4);

	Ch(1, 3) = Ch(0, 2) / R;
	Ch(3, 1) = Ch(1, 3);

	Ch(1, 1) = (PI*nt / a)*(R * R / (C3333*A) + a * a * c * c * (8 * PI * PI * nt * nt + 3) / (6 * C3333*I1)
		+ (c*c*c*c * (8 * PI * PI * nt * nt - 3) + 18 * R * R * (R * R - c * c)) / (6 * C3333*I2)
		+ (a * a * S1313*chi11 + c * c * S2323*chi22) / A
		+ R * R * c * c * (8 * PI * PI * nt * nt + 33) / (6 * J_th));

	Ch(1, 2) = -3 * R*Ch(1, 4) + 2 * PI*nt*R*c*a / (C3333*I1);
	Ch(2, 1) = Ch(1, 2);

	Ch(1, 4) = (PI*nt*c / (2 * a))*(-a * a / (C3333*I1) + (2 * R * R + c * c) / (C3333*I2) - R * R / J_th);
	Ch(4, 1) = Ch(1, 4);

	Ch(1, 5) = (2 * PI*nt*R / a)*((R * R - c * c) / (C3333*I2) + 2 * c * c / J_th);
	Ch(5, 1) = Ch(1, 5);

	Ch(2, 2) = (PI*nt / a)*(2 * c * c / (C3333*A) + a * a * R * R / (C3333*I1) 
		+ 3 * R * R * c * c / (C3333*I2) + 2 * R * R * S2323*chi22 / A + 3 * R * R * R * R / J_th);

	Ch(2, 4) = -R*Ch(1, 2) / (PI*nt*c * c) - 2 * PI*nt*R*a / (C3333*I1);
	Ch(4, 2) = Ch(2, 4);

	Ch(2, 5) = (2 * PI*nt*R * R * c / a)*(-1 / (C3333*I2) + 1 / J_th);
	Ch(5, 2) = Ch(2, 5);

	Ch(3, 3) = -Ch(2, 4) / R;

	Ch(4, 4) = -Ch(2, 4) / R;

	Ch(5, 5) = (2 * PI*nt / a)*(R * R / (C3333*I2) + c * c / J_th);

	// Stiffness matrix(single helix)
	Matrix Kh = invert6x6(Ch);

	// Stiffness matrix(assembly)
	D = nh*Kh;

	double deltas = sqrt((2 * PI*R)*(2 * PI*R) + p*p);
	double V = deltas*PI*r*r;
	double m = rho*V;
	double rhoL = m / p;
	double Ac = 4 * PI*R*r;
	double rhoE = m / (Ac*p);
	double Jp = rhoE*PI*0.5*(pow(R + r, 4) - pow(R - r, 4));
	
	Mr(0, 0) = rhoL;
	Mr(1, 1) = rhoL;
	Mr(2, 2) = rhoL;
	
	Jr(0, 0) = Jp / 2;
	Jr(1, 1) = Jp / 2;
	Jr(2, 2) = Jp;
}
