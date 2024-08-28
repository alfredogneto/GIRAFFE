class Matrix;

class SurfaceData
{
public:
	SurfaceData(int n_solutions);
	~SurfaceData();
	int n_sol;
	//Variaveis que descrevem as superficies envolvidas - cada uma delas e um vetor com 'n_solutions' posições, para descrever diversas posições de interesse da superficie em questão
	Matrix** G_p;
	Matrix** t1_p;
	Matrix** t2_p;
	Matrix** n_p;
	Matrix** G_i;
	Matrix** t1_i;
	Matrix** t2_i;
	Matrix** n_i;
	Matrix** G_ip;
};