#include "Environment.h"
#include <math.h>

#include "CoordinateSystem.h"

#include"Database.h"
#include "IO.h"
//Global variables
extern 
Database db;

#define PI 3.1415926535897932384626433832795

Environment::Environment(void)
{
	n_current_points = 0;
	n_wind_points = 0;
	current = NULL;
	wind = NULL;
	
	G = Matrix(3);
	surface_position = Matrix(3);
	reference_position = Matrix(3);
	transform3 = Matrix(3, 3);

	//Variables booleanas for controlar the that exists of data ambientais
	g_exist = false;
	ocean_data_exist = false;
	wind_data_exist = false;

	transform3_calculated = false;

	cs = 0;
}

Environment::~Environment(void)
{
	flush();
} 
//Sets the number of points of data of correnteza maritima
void Environment::SetNCurrentPoints(int value)
{
	flush();
	n_current_points = value;
	current = new double*[n_current_points];
	for (int i = 0; i < n_current_points; i++)
		current[i] = new double[3];
}
//Sets the number of points of data of the wind
void Environment::SetNWindPoints(int value)
{
	flush();
	n_wind_points = value;
	wind = new double*[n_wind_points];
	for (int i = 0; i < n_wind_points; i++)
		wind[i]= new double[8];
}
//Deletes alocações of the memory
void Environment::flush()
{
	if (current != NULL && n_current_points != 0)
	{
		for (int i = 0; i < n_current_points; i++)
				delete [] current[i];
			delete [] current;
	}
	current = NULL;
	if (wind != NULL && n_wind_points != 0)
	{
		for (int i = 0; i < n_wind_points; i++)
				delete[] wind[i];
		delete[] wind;
	}
	wind = NULL;
}
//Sets in the linha "line" the variables referring to the correnteza
void Environment::SetCurrentData(int line,double depth,double speed, double alpha)
{
	current[line][0] = depth;
	current[line][1] = speed;
	current[line][2] = alpha*PI/180;
}
//Sets in the linha "line" the variables referring to to the wind
void Environment::SetWindData(int line, double time, double wind_speed, double delta, double wind_vert,
	double HSHR, double VSHR, double Lin_VSHR, double gust_speed)
{
	wind[line][0] = time;
	wind[line][1] = wind_speed;
	wind[line][2] = delta*PI / 180;
	wind[line][3] = wind_vert;
	wind[line][4] = HSHR;
	wind[line][5] = VSHR;
	wind[line][6] = Lin_VSHR;
	wind[line][7] = gust_speed;
}
//Performs interpolation linear and returns the value of the velocity in certa profundidade
Matrix Environment::VelocityAt(double dep)
{
	Matrix ret_vel(3);
	int index = 0;
	//Check of the H - if is outside of the range of data of the corrente, returns a value null of velocity
	if (dep <= current[0][0])
	{
		ret_vel(0, 0) = current[0][1] * cos(current[0][2]);
		ret_vel(1, 0) = current[0][1] * sin(current[0][2]);
		return ret_vel;
	}
	else
	{
		if (dep >= current[n_current_points - 1][0])
		{
			ret_vel(0, 0) = current[n_current_points - 1][1] * cos(current[n_current_points - 1][2]);
			ret_vel(1, 0) = current[n_current_points - 1][1] * sin(current[n_current_points - 1][2]);
			return ret_vel;
		}
		//If not, does the interpolation
		else
		{
			while (dep > current[index][0])
				index++;
			if (index != 0)
			{
				//Interpolation linear
				//O point of interest will be between the indices index and index -1
				double factor = (dep - current[index - 1][0]) / (current[index][0] - current[index - 1][0]);
				double interp_speed = current[index - 1][1] + factor*(current[index][1] - current[index - 1][1]);
				double interp_alpha = current[index - 1][2] + factor*(current[index][2] - current[index - 1][2]);
				ret_vel(0, 0) = interp_speed*cos(interp_alpha);
				ret_vel(1, 0) = interp_speed*sin(interp_alpha);
			}
			else
			{
				db.myprintf("\nError in interpolation of sea current speed\n");
			}
			return ret_vel;
		}
	}
}
//Returns the Depth in the due index
double Environment::GetDepthAtIndex(int index)
{
	return current[index][0];
}
//Returns the Speed in the due index
double Environment::GetSpeedAtIndex(int index)
{
	return current[index][1];
}
//Returns the Angle in the due index
double Environment::GetAngleAtIndex(int index)
{
	return current[index][2];
}
//Performs calculation of the velocity horizontal
double Environment::VHor(double proj, double zlocal, double zhub, double wind_speed, double HSHR, double Lin_VSHR, double VSHR, double Vgust)
{
	double vret = wind_speed*(pow(((zlocal + zhub) / zhub), VSHR) + HSHR*proj + Lin_VSHR*zlocal) + Vgust;
	return vret;
	//Expressão of the aerodyn
	//V1 = V_tmp * ((InputPosition(3) / RefHt) ** VShr_tmp + (HShr_tmp * (InputPosition(2) * CosDelta + InputPosition(1) * SinDelta) + VLinShr_tmp * (InputPosition(3) - RefHt)) / RefWid) + VGUST_tmp
}

//Performs interpolation linear in the time
Matrix Environment::TimeAt(double time)
{
	Matrix ret_vel(8);
	int index = 0;
	//Check of the time - if is outside of the range of data of the wind, returns a value null of velocity
	if (time <= wind[0][0])
	{
		ret_vel(0, 0) = wind[0][0];
		ret_vel(1, 0) = wind[0][1];
		ret_vel(2, 0) = wind[0][2];
		ret_vel(3, 0) = wind[0][3];
		ret_vel(4, 0) = wind[0][4];
		ret_vel(5, 0) = wind[0][5];
		ret_vel(6, 0) = wind[0][6];
		ret_vel(7, 0) = wind[0][7];
		return ret_vel;
	}
	else
	{
		if (time >= wind[n_wind_points - 1][0])
		{
			ret_vel(0, 0) = wind[n_wind_points - 1][0];
			ret_vel(1, 0) = wind[n_wind_points - 1][1];
			ret_vel(2, 0) = wind[n_wind_points - 1][2];
			ret_vel(3, 0) = wind[n_wind_points - 1][3];
			ret_vel(4, 0) = wind[n_wind_points - 1][4];
			ret_vel(5, 0) = wind[n_wind_points - 1][5];
			ret_vel(6, 0) = wind[n_wind_points - 1][6];
			ret_vel(7, 0) = wind[n_wind_points - 1][7];
			return ret_vel;
		}
		//If not, does the interpolation
		else
		{
			while (time > wind[index][0])
				index++;
			if (index != 0)
			{
				//Interpolation linear
				//O point of interest will be between the indices index and index -1
				double factor = (time - wind[index - 1][0]) / (wind[index][0] - wind[index - 1][0]);
				double int_time = wind[index - 1][0] + factor*(wind[index][0] - wind[index - 1][0]);
				double int_windhor = wind[index - 1][1] + factor*(wind[index][1] - wind[index - 1][1]);
				double int_delta = wind[index - 1][2] + factor*(wind[index][2] - wind[index - 1][2]);
				double int_windvert = wind[index - 1][3] + factor*(wind[index][3] - wind[index - 1][3]);
				double int_HSHR = wind[index - 1][4] + factor*(wind[index][4] - wind[index - 1][4]);
				double int_VSHR = wind[index - 1][5] + factor*(wind[index][5] - wind[index - 1][5]);
				double int_VLIN = wind[index - 1][6] + factor*(wind[index][6] - wind[index - 1][6]);
				double int_VGUST = wind[index - 1][7] + factor*(wind[index][7] - wind[index - 1][7]);
				ret_vel(0, 0) = int_time;
				ret_vel(1, 0) = int_windhor;
				ret_vel(2, 0) = int_delta;
				ret_vel(3, 0) = int_windvert;
				ret_vel(4, 0) = int_HSHR;
				ret_vel(5, 0) = int_VSHR;
				ret_vel(6, 0) = int_VLIN;
				ret_vel(7, 0) = int_VGUST;
			}
			else
			{
				db.myprintf("\nError in interpolation of wind speed\n");
			}
			return ret_vel;
		}
	}
}
//Performs interpolation linear and returns the value of the velocity
Matrix Environment::WindVelocityAt(Matrix pos, double time)
{
	//Interpolation in the time
	Matrix vel_info(8);		//info already interpolada in the time
	Matrix local(3);		//position in the system local
	Matrix ref(3);			//reference in the system local
	Matrix V(3);			//vector velocity of the wind in the system local
	Matrix ret_vel(3);		//return of velocity in the system global
	vel_info = TimeAt(time);

	//Calculation of the transformation matrix of coordinates
	if (transform3_calculated == false)
	{
		//Conversion of the global coordinate system for the local
		Matrix e1(3);
		e1(0, 0) = 1.0;
		Matrix e2(3);
		e2(1, 0) = 1.0;
		Matrix e3(3);
		e3(2, 0) = 1.0;
		Matrix e1r = *db.CS[cs - 1]->E1;
		Matrix e2r = *db.CS[cs - 1]->E2;
		Matrix e3r = *db.CS[cs - 1]->E3;
		//Fills the transformation matrix of coordinates
		transform3(0, 0) = dot(e1r, e1);
		transform3(0, 1) = dot(e1r, e2);
		transform3(0, 2) = dot(e1r, e3);

		transform3(1, 0) = dot(e2r, e1);
		transform3(1, 1) = dot(e2r, e2);
		transform3(1, 2) = dot(e2r, e3);

		transform3(2, 0) = dot(e3r, e1);
		transform3(2, 1) = dot(e3r, e2);
		transform3(2, 2) = dot(e3r, e3);

		transform3_calculated = true;
	}

	//Conversion global-local
	local = transform3*pos;
	ref = transform3*reference_position;

	//Pegando data of the vel_info
	double wind_speed = vel_info(1, 0);
	double delta = vel_info(2, 0);
	double v_vertical = vel_info(3,0);
	double HSHR = vel_info(4, 0);
	double VSHR = vel_info(5, 0);
	double Lin_VSHR = vel_info(6, 0);
	double Vgust = vel_info(7, 0);
	
	//Translation - information of the reference
	local = local - ref;

	//Interpolation in the espaço
	//orientation of the plane of reference
	Matrix t(3);
	t(0, 0) = sin(delta);
	t(1, 0) = cos(delta);
	t(2, 0) = 0;
	//dot product local*t
	double proj = dot(t, local);
	double velh = VHor(proj, local(2, 0), ref(2, 0), wind_speed, HSHR, Lin_VSHR, VSHR, Vgust);

	// assembly of the ret_vel  
	ret_vel(0, 0) = velh * cos(delta);
	ret_vel(1, 0) = -velh * sin(delta);
	ret_vel(2, 0) = v_vertical;
	ret_vel = transp(transform3)*ret_vel;

	return ret_vel;
}
//Returns the Time in the due index
double Environment::GetTimeAt(int index)
{
	return wind[index][0];
}
//Returns the Wind Speed in the due index
double Environment::GetWindSpeedAt(int index)
{
	return wind[index][1];
}
//Returns the Delta in the due index
double Environment::GetDeltaAt(int index)
{
	return wind[index][2];
}
//Returns the Wind Vert in the due index
double Environment::GetWindVertAt(int index)
{
	return wind[index][3];
}
//Returns the Shear linear horizontal in the due index
double Environment::GetHorizShearAt(int index)
{
	return wind[index][4];
}
//Returns the Power Law Shear vertical in the due index
double Environment::GetVertPLShearAt(int index)
{
	return wind[index][5];
}
//Returns the Shear Linear vertical in the due index
double Environment::GetVertLinShearAt(int index)
{
	return wind[index][6];
}
//Returns the Gust Speed in the due index
double Environment::GetGustSpeedAt(int index)
{
	return wind[index][7];
}

bool Environment::Read(FILE *f)
{
	char s[1000];			//saves words-key lidas and values lidos
	fpos_t pos;				//variable that saves point of the stream of reading
	bool any_read = true;	//marca that some of the words-key was encontrada and realizada the reading

	while (any_read == true)
	{
		fgetpos(f, &pos);	//Saves the position (stream) of the start of the reading
		any_read = false;	//indicates that still not leu nada
		TryComment(f);
		//////////////////////////////////////////////////////
		//Reading of the GravityData
		fscanf(f, "%s", s);
		if (!strcmp(s, "GravityData") && any_read == false)
		{
			g_exist = true;
			fscanf(f, "%s", s);
			if (!strcmp(s, "G"))
			{
				fscanf(f, "%s", s);
				G(0, 0) = atof(s);
				fscanf(f, "%s", s);
				G(1, 0) = atof(s);
				fscanf(f, "%s", s);
				G(2, 0) = atof(s);

				//Reading of the BoolTable of the G
				fscanf(f, "%s", s);
				if (!strcmp(s, "BoolTable"))
				{
					bool_g.Read(f);
				}
				else
					return false;
			}
			any_read = true;		//Marca that was realizada reading
		}
		//////////////////////////////////////////////////////
		//Reading of the OceanData
		if (!strcmp(s, "OceanData") && any_read == false)
		{
			ocean_data_exist = true;
			//Reading of the RhoFluid
			fscanf(f, "%s", s);
			if (!strcmp(s, "RhoFluid"))
			{
				fscanf(f, "%s", s);
				rho_fluid = atof(s);
			}
			else
				return false;

			//Reading of the SurfacePosition	
			fscanf(f, "%s", s);
			if (!strcmp(s, "SurfacePosition"))
			{
				fscanf(f, "%s", s);
				surface_position(0, 0) = atof(s);
				fscanf(f, "%s", s);
				surface_position(1, 0) = atof(s);
				fscanf(f, "%s", s);
				surface_position(2, 0) = atof(s);
			}
			else
				return false;
			TryComment(f);
			//Reading of the corrente maritima
			fscanf(f, "%s", s);
			if (!strcmp(s, "SeaCurrent"))
			{
				//Reading of the number of points of the corrente (N)
				fscanf(f, "%s", s);
				if (!strcmp(s, "N"))
				{
					fscanf(f, "%s", s);
					n_current_points = atoi(s);
				}
				else
					return false;

				//Reading of the BoolTable
				fscanf(f, "%s", s);
				if (!strcmp(s, "BoolTable"))
				{
					bool_current.Read(f);
				}
				else
					return false;

				SetNCurrentPoints(n_current_points);
				//Reading of the points of the corrente
				for (int i = 0; i < n_current_points; i++)
				{
					TryComment(f);
					//Reading of each point of the corrente
					double tempD, tempSpeed, tempAngle;
					//Depth
					fscanf(f, "%s", s);
					if (!strcmp(s, "Depth"))
					{
						fscanf(f, "%s", s);
						tempD = atof(s);
					}
					else
						return false;

					//Speed
					fscanf(f, "%s", s);
					if (!strcmp(s, "Speed"))
					{
						fscanf(f, "%s", s);
						tempSpeed = atof(s);
					}
					else
						return false;

					//Angle
					fscanf(f, "%s", s);
					if (!strcmp(s, "Angle"))
					{
						fscanf(f, "%s", s);
						tempAngle = atof(s);
					}
					else
						return false;

					//Grava in the vector
					SetCurrentData(i, tempD, tempSpeed, tempAngle);
				}
			}
			else
				return false;
			any_read = true;		//Marca that was realizada reading
		}
		//////////////////////////////////////////////////////
		//Reading of the WindData
		if (!strcmp(s, "WindData") && any_read == false)
		{
			wind_data_exist = true;
			//Insert data of reading here
			//Reading of the Rho_air
			fscanf(f, "%s", s);
			if (!strcmp(s, "RhoAir"))
			{
				fscanf(f, "%s", s);
				rho_air = atof(s);
			}
			else
				return false;

			//Reading of the ReferencePosition	
			fscanf(f, "%s", s);
			if (!strcmp(s, "ReferencePosition"))
			{
				fscanf(f, "%s", s);
				reference_position(0, 0) = atof(s);
				fscanf(f, "%s", s);
				reference_position(1, 0) = atof(s);
				fscanf(f, "%s", s);
				reference_position(2, 0) = atof(s);
			}
			else
				return false;
			//Reading of the Coordinate system
			fscanf(f, "%s", s);
			if (!strcmp(s, "CS"))
			{
				fscanf(f, "%s", s);
				cs = atoi(s);
			}
			else
				return false;
			//Reading of the wind
			fscanf(f, "%s", s);
			if (!strcmp(s, "Wind"))
			{
				//Reading of the number of points of the wind (N)
				fscanf(f, "%s", s);
				if (!strcmp(s, "N"))
				{
					fscanf(f, "%s", s);
					n_wind_points = atoi(s);
				}
				else
					return false;
				SetNWindPoints(n_wind_points);
				//Reading of the points of the wind
				for (int i = 0; i < n_wind_points; i++)
				{
					//Reading of each point of the corrente
					double temp_time, temp_WindSpeed, temp_Delta, temp_WindVert, temp_HSHR, temp_VSHR, temp_VLIN, temp_VGUST;

					//Time
					fscanf(f, "%s", s);
					if (!strcmp(s, "Time"))
					{
						fscanf(f, "%s", s);
						temp_time = atof(s);
					}
					else
						return false;
					//Wind Speed
					fscanf(f, "%s", s);
					if (!strcmp(s, "WindSpeed"))
					{
						fscanf(f, "%s", s);
						temp_WindSpeed = atof(s);
					}
					else
						return false;
					//Delta
					fscanf(f, "%s", s);
					if (!strcmp(s, "Delta"))
					{
						fscanf(f, "%s", s);
						temp_Delta = atof(s);
					}
					else
						return false;
					//Wind Vert
					fscanf(f, "%s", s);
					if (!strcmp(s, "VerticalSpeed"))
					{
						fscanf(f, "%s", s);
						temp_WindVert = atof(s);
					}
					else
						return false;
					//Horizontal Shear
					fscanf(f, "%s", s);
					if (!strcmp(s, "HSHR"))
					{
						fscanf(f, "%s", s);
						temp_HSHR = atof(s);
					}
					else
						return false;
					//Vertical Shear
					fscanf(f, "%s", s);
					if (!strcmp(s, "VSHR"))
					{
						fscanf(f, "%s", s);
						temp_VSHR = atof(s);
					}
					else
						return false;
					//Linear Vertical Shear
					fscanf(f, "%s", s);
					if (!strcmp(s, "VLIN"))
					{
						fscanf(f, "%s", s);
						temp_VLIN = atof(s);
					}
					else
						return false;
					//Gust Speed
					fscanf(f, "%s", s);
					if (!strcmp(s, "VGUST"))
					{
						fscanf(f, "%s", s);
						temp_VGUST = atof(s);
					}
					else
						return false;
					//Check of the VSHR and VLIN
					if (temp_VSHR != 0.0 && temp_VLIN != 0.0)
					{
						db.myprintf("Error in WindData parameters: VSHR or VLIN has to be null!\n");
						return false;
					}
						
					//Grava in the vector
					SetWindData(i, temp_time, temp_WindSpeed, temp_Delta, temp_WindVert, temp_HSHR, temp_VSHR, temp_VLIN, temp_VGUST);
				}
			}
			else
				return false;
			any_read = true;		//Marca that was realizada reading
		}
		//////////////////////////////////////////////////////

		//If not there was reading of nada - sign of that already and a keyword not pertencente to the bloco "Environment". 
		//With this, saves the position of reading saves previous and for the reading.
		if (any_read == false)
			fsetpos(f, &pos);
	}//while
	return true;
}

void Environment::Write(FILE *f)
{
	if (g_exist == true)
	{
		fprintf(f, "GravityData\nG\t%.6e\t%.6e\t%.6e\t", G(0, 0), G(1, 0), G(2, 0));
		bool_g.Write(f);
	}
	if (ocean_data_exist == true)
	{
		fprintf(f, "OceanData\n");
		fprintf(f, "RhoFluid\t%.6e\tSurfacePosition\t%.6e\t%.6e\t%.6e\nSeaCurrent\tN\t%d\t",
			rho_fluid, surface_position(0, 0), surface_position(1, 0), surface_position(2, 0), n_current_points);
		bool_current.Write(f);
		for (int i = 0; i < n_current_points; i++)
		{
			fprintf(f, "Depth\t%lf\tSpeed\t%lf\tAngle\t%lf\n", GetDepthAtIndex(i), GetSpeedAtIndex(i), GetAngleAtIndex(i)*180.0 / PI);
		}
	}
	if (wind_data_exist == true)
	{
		fprintf(f, "WindData\n");
		fprintf(f, "RhoAir\t%.6e\tReferencePosition\t%.6e\t%.6e\t%.6e\tCS\t%d\nWind\tN\t%d\t",
			rho_air, reference_position(0, 0), reference_position(1, 0), reference_position(2, 0), cs, n_wind_points);
		
		for (int i = 0; i < n_wind_points; i++)
		{
			fprintf(f, "Time\t%lf\tWindSpeed\t%lf\tDelta\t%lf\tVerticalSpeed\t%lf\tHSHR\t%lf\tVSHR\t%lf\tVLIN\t%lf\tVGUST\t%lf\n", GetTimeAt(i), GetWindSpeedAt(i), GetDeltaAt(i)*180.0/PI, GetWindVertAt(i), GetHorizShearAt(i), GetVertPLShearAt(i), GetVertLinShearAt(i), GetGustSpeedAt(i));
		}
	}
	
}