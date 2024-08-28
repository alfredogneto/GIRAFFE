#include "Particle.h"
#include"Database.h"

//Variáveis globais
extern
Database db;

#define PI 3.1415926535897932384626433832795

Particle::Particle()
{
}

Particle::~Particle()
{
}

//Escreve no monitor da partícula
void Particle::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tkinetic\tang_mom_x\tang_mom_y\tang_mom_z\tang_mom_norm\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", time,
		kinetic_energy, angular_momentum[0], angular_momentum[1], angular_momentum[2], angular_momentum_mag);
}
void Particle::AllocMassFragmented(int nA, int nB)
{
	FreeMassFragmented();
	invMAA = new Matrix(nA, nA);
	MAB = new Matrix(nA, nB);
	MBB = new Matrix(nB, nB);
	GLA = new int[nA];
	for (int i = 0; i < nA; i++)
		GLA[i] = 0;
	GLB = new int[nB];
	for (int i = 0; i < nB; i++)
		GLB[i] = 0;
	allocedMassFragmented = true;
}
void Particle::FreeMassFragmented()
{
	if (allocedMassFragmented)
	{
		delete invMAA;
		delete MAB;
		delete MBB;
		delete[]GLA;
		delete[]GLB;
		allocedMassFragmented = false;
	}
}

