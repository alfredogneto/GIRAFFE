#pragma once
#include <stdio.h>
#include <vector>
using namespace std;

class Matrix;
class SurfacePairGeneralContact;

class ContactParticleBoundary
{
public:
	ContactParticleBoundary();
	virtual ~ContactParticleBoundary();

	int index1;				//Particle 1 - index
	int index2;				//Boundary 2 - index
	int sub_index1;			//Particle 1 - sub_index
	int sub_index2;			//Boundary 2 - sub_index
	bool cur_active;		//Flag active/unnactive - current status
	bool prev_active;		//Flag active/unnactive - previous status (of current time-step)
	 
	Matrix* I3;

	virtual void PreCalc() = 0;
	virtual void ProcessSurfacePairs() = 0;
	virtual void MountGlobal() = 0;
	virtual void MountGlobalExplicit() = 0;		
	virtual void ProcessContactHierarchy() = 0;
	//Explicit
	virtual void FinalProcessSurfacePairsExplicit(double t) = 0;

	void MountContacts();
	void MountContactsExplicit(double t);
	void SaveConfiguration();
	bool HaveErrors();
	void Clear();
	bool NightOwlContact();
	void WriteVTK_XMLForces(FILE *f);
	double TimeStepControl(double kin);
	double EvaluateContactEnergy();

	//Vector of contact pairs (to be alloced according to the particle - particle type)
	vector<SurfacePairGeneralContact*> contact_pairs;
};

