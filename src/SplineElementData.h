class Matrix;

class SplineElementData
{
public:
	SplineElementData(int n_solutions);
	~SplineElementData();
	int n_sol;
	//Variaveis que descrevem as superficies envolvidas - cada uma delas e um vetor com 'n_solutions' posi��es, para descrever diversas posi��es de interesse da superficie em quest�o
	Matrix** G_p;

};