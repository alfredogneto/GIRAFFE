#include "Environment.h"
#include <math.h>

#include "CoordinateSystem.h"

#include"Database.h"
#include "IO.h"
//Variaveis globais
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

	//Variaveis booleanas para controlar o que existe de dados ambientais
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
//Seta o numero de pontos de dados de correnteza maritima
void Environment::SetNCurrentPoints(int value)
{
	flush();
	n_current_points = value;
	current = new double*[n_current_points];
	for (int i = 0; i < n_current_points; i++)
		current[i] = new double[3];
}
//Seta o numero de pontos de dados do vento
void Environment::SetNWindPoints(int value)
{
	flush();
	n_wind_points = value;
	wind = new double*[n_wind_points];
	for (int i = 0; i < n_wind_points; i++)
		wind[i]= new double[8];
}
//Deleta alocações da memória
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
//Seta na linha "line" as variaveis referentes a correnteza
void Environment::SetCurrentData(int line,double depth,double speed, double alpha)
{
	current[line][0] = depth;
	current[line][1] = speed;
	current[line][2] = alpha*PI/180;
}
//Seta na linha "line" as variaveis referentes ao vento
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
//Realiza interpolação linear e retorna o valor da velocidade em certa profundidade
Matrix Environment::VelocityAt(double dep)
{
	Matrix ret_vel(3);
	int index = 0;
	//Verificação do H - se estiver fora do range de dados da corrente, retorna um valor nulo de velocidade
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
		//Se não, faz a interpolação
		else
		{
			while (dep > current[index][0])
				index++;
			if (index != 0)
			{
				//Interpolação linear
				//O ponto de interesse estara entre os indices index e index -1
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
//Retorna o Depth no devido index
double Environment::GetDepthAtIndex(int index)
{
	return current[index][0];
}
//Retorna o Speed no devido index
double Environment::GetSpeedAtIndex(int index)
{
	return current[index][1];
}
//Retorna o Angle no devido index
double Environment::GetAngleAtIndex(int index)
{
	return current[index][2];
}
//Realiza calculo da velocidade horizontal
double Environment::VHor(double proj, double zlocal, double zhub, double wind_speed, double HSHR, double Lin_VSHR, double VSHR, double Vgust)
{
	double vret = wind_speed*(pow(((zlocal + zhub) / zhub), VSHR) + HSHR*proj + Lin_VSHR*zlocal) + Vgust;
	return vret;
	//Expressão do aerodyn
	//V1 = V_tmp * ((InputPosition(3) / RefHt) ** VShr_tmp + (HShr_tmp * (InputPosition(2) * CosDelta + InputPosition(1) * SinDelta) + VLinShr_tmp * (InputPosition(3) - RefHt)) / RefWid) + VGUST_tmp
}

//Realiza interpolação linear no tempo
Matrix Environment::TimeAt(double time)
{
	Matrix ret_vel(8);
	int index = 0;
	//Verificação do tempo - se estiver fora do range de dados do vento, retorna um valor nulo de velocidade
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
		//Se não, faz a interpolação
		else
		{
			while (time > wind[index][0])
				index++;
			if (index != 0)
			{
				//Interpolação linear
				//O ponto de interesse estara entre os indices index e index -1
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
//Realiza interpolação linear e retorna o valor da velocidade
Matrix Environment::WindVelocityAt(Matrix pos, double time)
{
	//Interpolação no tempo
	Matrix vel_info(8);		//info ja interpolada no tempo
	Matrix local(3);		//posição no sistema local
	Matrix ref(3);			//referência no sistema local
	Matrix V(3);			//vetor velocidade do vento no sistema local
	Matrix ret_vel(3);		//retorno de velocidade no sistema global
	vel_info = TimeAt(time);

	//Calculo da matriz de transformação de coordenadas
	if (transform3_calculated == false)
	{
		//Conversão do sistema de coordenadas global para o local
		Matrix e1(3);
		e1(0, 0) = 1.0;
		Matrix e2(3);
		e2(1, 0) = 1.0;
		Matrix e3(3);
		e3(2, 0) = 1.0;
		Matrix e1r = *db.CS[cs - 1]->E1;
		Matrix e2r = *db.CS[cs - 1]->E2;
		Matrix e3r = *db.CS[cs - 1]->E3;
		//Preenche a matriz de transformação de coordenadas
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

	//Conversão global-local
	local = transform3*pos;
	ref = transform3*reference_position;

	//Pegando dados do vel_info
	double wind_speed = vel_info(1, 0);
	double delta = vel_info(2, 0);
	double v_vertical = vel_info(3,0);
	double HSHR = vel_info(4, 0);
	double VSHR = vel_info(5, 0);
	double Lin_VSHR = vel_info(6, 0);
	double Vgust = vel_info(7, 0);
	
	//Translação - informação da referência
	local = local - ref;

	//Interpolação no espaço
	//orientação do plano de referência
	Matrix t(3);
	t(0, 0) = sin(delta);
	t(1, 0) = cos(delta);
	t(2, 0) = 0;
	//produto escalar local*t
	double proj = dot(t, local);
	double velh = VHor(proj, local(2, 0), ref(2, 0), wind_speed, HSHR, Lin_VSHR, VSHR, Vgust);

	// montagem do ret_vel  
	ret_vel(0, 0) = velh * cos(delta);
	ret_vel(1, 0) = -velh * sin(delta);
	ret_vel(2, 0) = v_vertical;
	ret_vel = transp(transform3)*ret_vel;

	return ret_vel;
}
//Retorna o Time no devido index
double Environment::GetTimeAt(int index)
{
	return wind[index][0];
}
//Retorna o Wind Speed no devido index
double Environment::GetWindSpeedAt(int index)
{
	return wind[index][1];
}
//Retorna o Delta no devido index
double Environment::GetDeltaAt(int index)
{
	return wind[index][2];
}
//Retorna o Wind Vert no devido index
double Environment::GetWindVertAt(int index)
{
	return wind[index][3];
}
//Retorna o Shear linear horizontal no devido index
double Environment::GetHorizShearAt(int index)
{
	return wind[index][4];
}
//Retorna o Power Law Shear vertical no devido index
double Environment::GetVertPLShearAt(int index)
{
	return wind[index][5];
}
//Retorna o Shear Linear vertical no devido index
double Environment::GetVertLinShearAt(int index)
{
	return wind[index][6];
}
//Retorna o Gust Speed no devido index
double Environment::GetGustSpeedAt(int index)
{
	return wind[index][7];
}

bool Environment::Read(FILE *f)
{
	char s[1000];			//salva palavras-chave lidas e valores lidos
	fpos_t pos;				//variavel que salva ponto do stream de leitura
	bool any_read = true;	//marca que alguma das palavras-chave foi encontrada e realizada a leitura

	while (any_read == true)
	{
		fgetpos(f, &pos);	//Salva a posição (stream) do inicio da leitura
		any_read = false;	//indica que ainda não leu nada
		TryComment(f);
		//////////////////////////////////////////////////////
		//Leitura do GravityData
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

				//Leitura do BoolTable do G
				fscanf(f, "%s", s);
				if (!strcmp(s, "BoolTable"))
				{
					bool_g.Read(f);
				}
				else
					return false;
			}
			any_read = true;		//Marca que foi realizada leitura
		}
		//////////////////////////////////////////////////////
		//Leitura do OceanData
		if (!strcmp(s, "OceanData") && any_read == false)
		{
			ocean_data_exist = true;
			//Leitura do RhoFluid
			fscanf(f, "%s", s);
			if (!strcmp(s, "RhoFluid"))
			{
				fscanf(f, "%s", s);
				rho_fluid = atof(s);
			}
			else
				return false;

			//Leitura do SurfacePosition	
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
			//Leitura da corrente maritima
			fscanf(f, "%s", s);
			if (!strcmp(s, "SeaCurrent"))
			{
				//Leitura do numero de pontos da corrente (N)
				fscanf(f, "%s", s);
				if (!strcmp(s, "N"))
				{
					fscanf(f, "%s", s);
					n_current_points = atoi(s);
				}
				else
					return false;

				//Leitura do BoolTable
				fscanf(f, "%s", s);
				if (!strcmp(s, "BoolTable"))
				{
					bool_current.Read(f);
				}
				else
					return false;

				SetNCurrentPoints(n_current_points);
				//Leitura dos pontos da corrente
				for (int i = 0; i < n_current_points; i++)
				{
					TryComment(f);
					//Leitura de cada ponto da corrente
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

					//Grava no vetor
					SetCurrentData(i, tempD, tempSpeed, tempAngle);
				}
			}
			else
				return false;
			any_read = true;		//Marca que foi realizada leitura
		}
		//////////////////////////////////////////////////////
		//Leitura do WindData
		if (!strcmp(s, "WindData") && any_read == false)
		{
			wind_data_exist = true;
			//Inserir dados de leitura aqui
			//Leitura do Rho_air
			fscanf(f, "%s", s);
			if (!strcmp(s, "RhoAir"))
			{
				fscanf(f, "%s", s);
				rho_air = atof(s);
			}
			else
				return false;

			//Leitura do ReferencePosition	
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
			//Leitura do Sistema de Coordenadas
			fscanf(f, "%s", s);
			if (!strcmp(s, "CS"))
			{
				fscanf(f, "%s", s);
				cs = atoi(s);
			}
			else
				return false;
			//Leitura do vento
			fscanf(f, "%s", s);
			if (!strcmp(s, "Wind"))
			{
				//Leitura do numero de pontos do vento (N)
				fscanf(f, "%s", s);
				if (!strcmp(s, "N"))
				{
					fscanf(f, "%s", s);
					n_wind_points = atoi(s);
				}
				else
					return false;
				SetNWindPoints(n_wind_points);
				//Leitura dos pontos do vento
				for (int i = 0; i < n_wind_points; i++)
				{
					//Leitura de cada ponto da corrente
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
					//Verificação do VSHR e VLIN
					if (temp_VSHR != 0.0 && temp_VLIN != 0.0)
					{
						db.myprintf("Error in WindData parameters: VSHR or VLIN has to be null!\n");
						return false;
					}
						
					//Grava no vetor
					SetWindData(i, temp_time, temp_WindSpeed, temp_Delta, temp_WindVert, temp_HSHR, temp_VSHR, temp_VLIN, temp_VGUST);
				}
			}
			else
				return false;
			any_read = true;		//Marca que foi realizada leitura
		}
		//////////////////////////////////////////////////////

		//Se não houve leitura de nada - sinal de que ja e uma palavra-chave não pertencente ao bloco "Environment". 
		//Com isso, salva a posição de leitura salva anterior e para a leitura.
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