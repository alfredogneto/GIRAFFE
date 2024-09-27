class Matrix;

class SplineElementData
{
public:
	SplineElementData(int n_solutions);
	~SplineElementData();
	int n_sol;
	//Variaveis que descrevem as superficies envolvidas - cada uma delas e um vetor com 'n_solutions' posições, para descrever diversas posições de interesse da superficie em questão
	Matrix** G_p;

};