#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class PostFiles
{
public:
	PostFiles();
	~PostFiles();

	double mag_factor;		//Fator de magnificação de deslocamentos

	int alloced_files;
	
	FILE **f_mesh_pvd;
	FILE **f_rendermesh_pvd;
	FILE **f_contactsurfaces_pvd;
	FILE **f_symbols_pvd;
	FILE **f_contactforces_pvd;
	FILE **f_rb_particles_pvd;
	FILE **f_constraints_pvd;
	FILE **f_forces_pvd;
	FILE **f_spline_pvd;
	FILE **f_renderspline_pvd;

	bool WriteMesh_flag;
	bool WriteRenderMesh_flag;
	bool WriteRigidContactSurfaces_flag;
	bool WriteFlexibleContactSurfaces_flag;
	bool WriteConstraints_flag;
	bool WriteForces_flag;
	bool WriteSpecialConstraints_flag;
	bool WriteContactForces_flag;
	bool WriteRenderRigidBodies_flag;
	bool WriteRenderParticles_flag;
	bool WriteSpline_flag;
	bool WriteRenderSpline_flag;
	
	bool Read(FILE *f);
	void Write(FILE *f);
	void AllocFiles(int n_files);
	void FlushFiles();
	void StartPostFiles(int sol_index);
	void EndPostFiles(int sol_index);
	void StartSinglePartFiles(int sol_index);
	void StartMultiplePartFiles(int sol_index);
	void UpdateSinglePartPostFiles(int sol_index, double time, int index);
	void UpdateMultiplePartPostFiles(int sol_index, double time, int index, int part);
	void EndSinglePartPostFiles(int sol_index);
	void EndMultiplePartPostFiles(int sol_index);
	
	void WriteVTK_XMLRender(int sol_index, double time, int index, int part);	//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLBase(int sol_index, double time, int index, int part);		//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLContact(int sol_index, double time, int index);			//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLSymbols(int sol_index, double time, int index);			//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLContactForces(int sol_index, double time, int index);		//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLRBParticles(int sol_index, double time, int index);		//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLConstraints(int sol_index, double time, int index);		//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLForces(int sol_index, double time, int index);				//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLSpline(int sol_index, double time, int index);				//Gerando arquivo para o Paraview (VTK) - formato XML
	void WriteVTK_XMLRenderSpline(int sol_index, double time, int index);		//Gerando arquivo para o Paraview (VTK) - formato XML

	void WriteConfigurationResults(int sol_index, double time, int index);

	void WriteSingleVertexPart(FILE *f);
};

