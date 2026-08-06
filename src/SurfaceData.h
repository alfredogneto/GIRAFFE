class Matrix;

class SurfaceData
{
public:
	SurfaceData(int n_solutions);
	~SurfaceData();
	int n_sol;
	//Variables that descrevem the surfaces envolvidas - each a delas and a vector with 'n_solutions' positions, for descrever diversas positions of interest of the surface in question
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