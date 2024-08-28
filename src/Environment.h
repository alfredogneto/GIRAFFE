#pragma once
#include <stdio.h>

#include "Matrix.h"
#include "BoolTable.h"

class Environment
{
public:
	Environment(void);
	~Environment(void);
	void flush();															//Deleta alocações da memória
	bool Read(FILE *f);
	void Write(FILE *f);

	//Variáveis booleanas para controlar o que existe de dados ambientais
	bool g_exist;
	bool ocean_data_exist;
	bool wind_data_exist;

	//BoolTables
	BoolTable bool_g;																//Sequencial de passos para o peso próprio
	BoolTable bool_current;															//Sequencial de passos para a corrente
	
	//Gravity data
	Matrix G;																//Campo gravitacional
	
	//Sea current data
	void SetNCurrentPoints(int value);										//Seta o número de pontos de dados de correnteza marítima
	void SetCurrentData(int line, double depth, double speed, double alpha);	//Seta na linha "line" as variáveis referente à correnteza
	Matrix VelocityAt(double H);											//Realiza interpolação linear e retorna o valor da velocidade em certa profundidade
	double GetDepthAtIndex(int index);										//Retorna o Depth no devido index								
	double GetSpeedAtIndex(int index);										//Retorna o Speed no devido index
	double GetAngleAtIndex(int index);										//Retorna o Angle no devido index
	double rho_fluid;														//Massa específica da água
	Matrix surface_position;												//Posição da lâmina d'água no ambiente
	int n_current_points;													//Número de pontos de dados da correnteza marítima
	double** current;														//Matriz que guarda os dados referentes à correnteza marítima
	
	//WindData
	Matrix TimeAt(double time);												//Realiza interpolação linear no tempo - retorna a matriz com todos os dados já interpolados no tempo
	Matrix WindVelocityAt(Matrix pos, double time);							//Realiza interpolação linear no espaço e tempo - retorna vetor velocidade no ponto de interesse
	double GetTimeAt(int index);											//Retorna o Time no devido index
	double GetWindSpeedAt(int index);										//Retorna o Wind Speed no devido index
	double GetDeltaAt(int index);											//Retorna o Delta no devido index
	double GetWindVertAt(int index);										//Retorna o Wind Vert no devido index
	double GetHorizShearAt(int index);										//Retorna o Shear Linear horizontal no devido index
	double GetVertPLShearAt(int index);										//Retorna o Power Law Shear vertical no devido index
	double GetVertLinShearAt(int index);									//Retorna o Shear Linear vertical no devido index
	double GetGustSpeedAt(int index);										//Retorna o Gust Speed no devido index

	void SetWindData(int line, double time, double wind_speed, double delta, double wind_vert, double HSHR, double VSHR, double Lin_VSHR, double gust_speed);		//Seta na linha "line" as variáveis referentes ao vento

	/*
	wind[line][0] = time;
	wind[line][1] = wind_speed;
	wind[line][2] = delta*PI / 180;
	wind[line][3] = wind_vert;
	wind[line][4] = HSHR;
	wind[line][5] = VSHR;
	wind[line][6] = Lin_VSHR;
	wind[line][7] = gust_speed;
	
	*/
	int cs;																	//ID do sistema de coordenadas para definir o vento
	
	double VHor(double proj, double zlocal, double zhub, double wind_speed, double HSHR, double Lin_VSHR, double VSHR, double Vgust);	//Realiza cálculo da velocidade horizontal
	double rho_air;															//Massa específica do ar
	Matrix reference_position;												//Posição da incidência do vento - sistema global
	Matrix transform3;														//Matriz de transformação de coordenadas 3x3 - local-global
	bool transform3_calculated;												//Variável booleana que indica que a matriz de transformação de coordenadas já foi calculada
	
	double** wind;															//Matriz que guarda os dados referentes ao vento (todos os dados de leitura - para cada instante de tempo)
	
	int n_wind_points;														//Número de pontos de dados do vento (número de instantes de tempo)
	void SetNWindPoints(int value);											//Seta o número de pontos de dados do ar (número de instantes - alocação de memória)
	
};
