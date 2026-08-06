class Matrix;

class SplineElementData
{
public:
	SplineElementData(int n_solutions);
	~SplineElementData();
	int n_sol;
	//Variables that descrevem the surfaces envolvidas - each a delas and a vector with 'n_solutions' positions, for descrever diversas positions of interest of the surface in question
	Matrix** G_p;

};