#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "BoundingVolume.h"
#define PI 3.1415926535897932384626433832795


class Geometry
{
public:
	Geometry();
	virtual ~Geometry();
	void Alloc();
	void Free();
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual void AllocSpecific() = 0;
	virtual void FreeSpecific() = 0;
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual bool Check() = 0;
	virtual void PreCalc() = 0;

	virtual void UpdateVariables() = 0;
	virtual void UpdateBoundingVolumes() = 0;				//Updates bounding volumes
	virtual void SaveLagrange() = 0;						//Salva variáveis

	virtual void SurfacePosition(Matrix* Gamma, double z, double th, bool next) = 0;	//returns the surface position on the vector Gamma, evaluated at (z,th) on current (i) or next (i+1) configuration
	virtual void NormalExt(double z, double th, Matrix* n, bool next) = 0;							//returns the exterior normal of the surface
	//Degeneration (true or false for first and second convective coordinates - when applicable for surfaces)
	bool deg_u1, deg_u2, deg_u1_u2;
	//Degenerated coordinates (only considered if boolean degeneration values are TRUE)
	double deg_u1_value, deg_u2_value, deg_u1_u2_relation;

	BoundingVolume* bv;		//Bounding volume
	float inc_len_factor;	//increase factor do BV
	float bv_offset;		//offset do BV
	float max[3], min[3];	//coordinates max and min - global CS - to guide the BodyGeometry Axis Aligned BV
	float largest_gnb;		//para auxiliar no cálculo do 'size' do bounding volume - Verlet/LinkedCells schemes

	int number;				//ID da Geometry
	int material;			//ID do material
	int super_node;			//Super node (quando aplicável)
	int mother_entity;		//ID da superfície mãe usada para degeneração (para casos de hierarquia de pares de contato)
	char* type_name;		//Nome do tipo de Geometry
	int *nodes;				//Lista de nós
	int n_nodes;			//Número de nós
	int **DOFs;				//Indica para a indexação de cada grau de liberdade, 1 ou 0, ativo ou inativo
	int nDOFs;				//Número do GL (global)
	int **GLs;				//Ponteiro para os GL globais utilizados

	double* body_mass;		//Mass of the associated body (pointer)

	//Variables to store DOFs and derivatives
	Matrix* d;
	Matrix* dui;
	Matrix* ddui;
};