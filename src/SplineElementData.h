#include "Matrix.h"
class SplineElementData
{
public:
	SplineElementData(int n_solutions);
	~SplineElementData();
	int n_sol;
	//Vari�veis que descrevem as superf�cies envolvidas - cada uma delas � um vetor com 'n_solutions' posi��es, para descrever diversas posi��es de interesse da superf�cie em quest�o
	Matrix** G_p;

};