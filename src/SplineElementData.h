#include "Matrix.h"
class SplineElementData
{
public:
	SplineElementData(int n_solutions);
	~SplineElementData();
	int n_sol;
	//Variáveis que descrevem as superfícies envolvidas - cada uma delas é um vetor com 'n_solutions' posições, para descrever diversas posições de interesse da superfície em questão
	Matrix** G_p;

};