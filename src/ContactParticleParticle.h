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
	bool *contact_detection; // Para cada iteração, essa variável guarda a informação se o contato entre as partículas já foi detectado (sem degeneração)
	bool *contact_detectionA; // Para cada iteração, essa variável guarda a informação se o contato entre as partículas já foi detectado projetando-se o nó da partícula A na partícula B
	bool *contact_detectionB; // Para cada iteração, essa variável guarda a informação se o contato entre as partículas já foi detectado projetando-se o nó da partícula B na partícula A
	bool *contact_detectionAB; // Para cada iteração, essa variável guarda a informação se o contato entre as partículas já foi detectado através da distância entre os nós das partículas A e B
};
