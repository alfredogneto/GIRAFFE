#pragma once
#include <vector>

using namespace std;

class Matrix;
class SurfacePairGeneralContact;

class ContactParticleParticle
{
public:
	ContactParticleParticle();
	virtual ~ContactParticleParticle();

	int index1;				//Particle 1 - index
	int index2;				//Particle 2 - index
	int sub_index1;			//Particle 1 - sub_index
	int sub_index2;			//Particle 2 - sub_index
	bool cur_active;		//Flag active/unnactive - current status
	bool prev_active;		//Flag active/unnactive - previous status (of current time-step)
	
	Matrix* I3;

	virtual void PreCalc() = 0;
	virtual void ProcessSurfacePairs() = 0;
	//Explicit
	virtual void FinalProcessSurfacePairsExplicit(double t) = 0;
	virtual void MountGlobal() = 0;
	virtual void MountGlobalExplicit() = 0;
	virtual void ProcessContactHierarchy() = 0;
	void MountContacts();
	void MountContactsExplicit(double t);
	void SaveConfiguration();
	bool HaveErrors();
	void Clear();
	bool NightOwlContact();
	void WriteVTK_XMLForces(FILE *f);
	double TimeStepControl(double kin);

	//Vector of contact pairs (to be alloced according to the particle - particle type)
	vector<SurfacePairGeneralContact*> contact_pairs;
	//Marina
	bool *contact_detection; // Para cada itera��o, essa vari�vel guarda a informa��o se o contato entre as part�culas j� foi detectado (sem degenera��o)
	bool *contact_detectionA; // Para cada itera��o, essa vari�vel guarda a informa��o se o contato entre as part�culas j� foi detectado projetando-se o n� da part�cula A na part�cula B
	bool *contact_detectionB; // Para cada itera��o, essa vari�vel guarda a informa��o se o contato entre as part�culas j� foi detectado projetando-se o n� da part�cula B na part�cula A
	bool *contact_detectionAB; // Para cada itera��o, essa vari�vel guarda a informa��o se o contato entre as part�culas j� foi detectado atrav�s da dist�ncia entre os n�s das part�culas A e B
};
