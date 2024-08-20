#pragma once
#include "BoundingVolume.h"
#include <vector>
#include <string>
#include <ctype.h>
#include <cstdio>
#include <array>
#include <omp.h>

using namespace std;

class Verlet
{
public:
	Verlet();
	~Verlet();
	bool Read(FILE *f);					//Leitura
	void Write(FILE *f);				//Gravação

	void AllocTables();
	void FreeTables();
	void RefreshVerletTables();
	void SetSamplingRate();

	//Verlet Tables for particles and sub-particles PP
	vector<std::array<int, 3>>* neighborhood_PP;
	//Verlet Tables for particles and sub-particles PB
	vector<std::array<int, 3>>* neighborhood_PB;
	//Verlet Tables for bodies and sub-bodies BOBO
	vector<std::array<int, 3>>* neighborhood_BOBO;
	//Verlet Tables for bodies and sub-bodies BOBO
	vector<std::array<int, 3>>* neighborhood_PBO;

	void ComputeVerletUpdateItems(float* C_i, float* C_j, float size_i, float size_j, float* l_i, float* l_j, int* Nsampling, bool* bool_near);


	float f1;	//cut-off maximum factor (for disregarding contact search)
	float f2;	//cut-off minimum factor (for re-sampling)
	//f1 > f2, f1 > 0 and f2 > 0

	//Combine with LinkedCells
	bool combine_linked_cells;
	//Run-time variables for sampling
	int acc_steps;
	int n_sampling;
	int MAX_SAMPLE_RATE;
};

