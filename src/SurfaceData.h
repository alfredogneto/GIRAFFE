class Matrix;

class SurfaceData
{
public:
	SurfaceData(int n_solutions);
	~SurfaceData();
	int n_sol;
	//Vari�veis que descrevem as superf�cies envolvidas - cada uma delas � um vetor com 'n_solutions' posi��es, para descrever diversas posi��es de interesse da superf�cie em quest�o
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