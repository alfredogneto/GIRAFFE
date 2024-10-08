#pragma once
class SSContactData;
class TimeStepControlData;
class Matrix;

class SurfacePairGeneralContact
{
public:
	SurfacePairGeneralContact();
	virtual ~SurfacePairGeneralContact();
	int index1;				//index 1
	int index2;				//index 2
	int sub_index1;			//sub index 1
	int sub_index2;			//sub index 2

	//Degenerations to employ while considering contact
	int deg_pointA;	//ID of the point in surface A to be considered in the degeneration (zero if no degeneration)
	int deg_curveA; //ID of the curve in surface A to be considered in the degeneration (zero if no degeneration)
	int deg_pointB; //ID of the point in surface B to be considered in the degeneration (zero if no degeneration)
	int deg_curveB; //ID of the curve in surface B to be considered in the degeneration (zero if no degeneration)

	//Variables for contact evaluation:
	SSContactData* cd;						//information of the contact between surfaces
	double* Rc;								//Vetor de esfor�os internos
	double** Kc;							//Matriz de rigidez

	//Time step control data
	TimeStepControlData* td;
	
	bool GetActive();
	void EvaluateInvertedHessian();			//Calcula a inversa da Hessiana
	void SetActive();
	void SetUnnactive();
	void Alloc();
	void Free();
	
	virtual void EvaluateNormalGap() = 0;
	virtual bool HaveErrors() = 0;
	virtual void SolveLCP() = 0;
	virtual void AllocSpecific() = 0;
	virtual void FreeSpecific() = 0;
	virtual void MountLocalContributions() = 0;

	virtual void MountLocalContributionsExplicit(double t) = 0;
	virtual void SetVariablesExplicit(double t) = 0;
	virtual void FinalUpdateExplicit(double t) = 0;

	virtual void PreCalc() = 0;
	virtual void SetVariables() = 0;					//Sets variables for AceGen codes interfaces
	virtual void HessianPhase1(Matrix& mHes) = 0;
	virtual void SurfacePoints() = 0;					//Sets GammaA and GammaB with contact positions on both surfaces
	virtual void Report() = 0;
	virtual void CompactReport() = 0;
	virtual void PredictorTimeStep(double kin) = 0;

	Matrix* GammaA;
	Matrix* GammaB;
	//Common AceGenPointers (to all surface pairs)
	double* fn;
	double* ft;
	double* vnrel;
	double* v;

	bool eligible;			//Set by the function that creates surface pairs - may change in each iteration of NR
	bool prev_eligible;		//Previous converged eligible


	bool *impactcontrol;
	
protected:
	bool active;			//Possible active/unnactive (due to neighboring contact)
};