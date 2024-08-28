#include "PostFiles.h"
#include"Database.h"
#include <direct.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Static.h"
#include "Dynamic.h"
#include "ExplicitDynamic.h"
#include "Modal.h"
#include "Element.h"
#include "RigidBody_1.h"
#include "Particle.h"
#include "Contact.h"
#include "Boundary.h"
#include "Geometry.h"
#include "GeneralContactSearch.h"
#include "Node.h"
#include "Encoding.h"
#include "Sphere.h"
#include "Spline.h"
#include "SpecialConstraint.h"
#include "Constraint.h"
#include "Load.h"

//Variáveis globais
extern
Database db;

PostFiles::PostFiles()
{
	mag_factor = 1;		//Fator de magnificação de deslocamentos
	alloced_files = 0;

	WriteMesh_flag = true;						//default
	WriteRenderMesh_flag = false;				
	WriteRigidContactSurfaces_flag = false;
	WriteFlexibleContactSurfaces_flag = false;
	WriteConstraints_flag = false;
	WriteForces_flag = false;
	WriteSpecialConstraints_flag = false;
	WriteContactForces_flag = false;
	WriteRenderRigidBodies_flag = false;
	WriteRenderParticles_flag = false;
	WriteSpline_flag = false;
	WriteRenderSpline_flag = false;

	f_mesh_pvd = NULL;
	f_rendermesh_pvd = NULL;
	f_contactsurfaces_pvd = NULL;
	f_symbols_pvd = NULL;
	f_contactforces_pvd = NULL;
	f_rb_particles_pvd = NULL;
	f_constraints_pvd = NULL;
	f_forces_pvd = NULL;
	f_spline_pvd = NULL;
	f_renderspline_pvd = NULL;
}

PostFiles::~PostFiles()
{
	FlushFiles();
}

void PostFiles::AllocFiles(int n_files)
{
	f_mesh_pvd = new FILE*[n_files];
	f_rendermesh_pvd = new FILE*[n_files];
	f_contactsurfaces_pvd = new FILE*[n_files];
	f_symbols_pvd = new FILE*[n_files];
	f_contactforces_pvd = new FILE*[n_files];
	f_rb_particles_pvd = new FILE*[n_files];
	f_constraints_pvd = new FILE*[n_files];
	f_forces_pvd = new FILE*[n_files];
	f_spline_pvd = new FILE*[n_files];
	f_renderspline_pvd = new FILE*[n_files];
	alloced_files = n_files;
}
void PostFiles::FlushFiles()
{
	if (alloced_files != 0)
	{
		if (f_mesh_pvd != NULL)
			delete[]f_mesh_pvd;
		if (f_rendermesh_pvd != NULL)
			delete[]f_rendermesh_pvd;
		if (f_contactsurfaces_pvd != NULL)
			delete[]f_contactsurfaces_pvd;
		if (f_symbols_pvd != NULL)
			delete[]f_symbols_pvd;
		if (f_contactforces_pvd != NULL)
			delete[]f_contactforces_pvd;
		if (f_rb_particles_pvd != NULL)
			delete[]f_rb_particles_pvd;
		if (f_constraints_pvd != NULL)
			delete[]f_constraints_pvd;
		if (f_forces_pvd != NULL)
			delete[]f_forces_pvd;
		if (f_spline_pvd != NULL)
			delete[]f_spline_pvd;
		if (f_renderspline_pvd != NULL)
			delete[]f_renderspline_pvd;

		f_mesh_pvd = NULL;
		f_rendermesh_pvd = NULL;
		f_contactsurfaces_pvd = NULL;
		f_symbols_pvd = NULL;
		f_contactforces_pvd = NULL;
		f_rb_particles_pvd = NULL;
		f_constraints_pvd = NULL;
		f_forces_pvd = NULL;
		f_spline_pvd = NULL;
		f_renderspline_pvd = NULL;
		alloced_files = 0;
	}
}

bool PostFiles::Read(FILE *f)
{
	char s[1000];
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "MagFactor"))
	{
		fscanf(f, "%s", s);
		mag_factor = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteMesh"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteMesh_flag = true;
		else
			WriteMesh_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteRenderMesh"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteRenderMesh_flag = true;
		else
			WriteRenderMesh_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteRigidContactSurfaces"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteRigidContactSurfaces_flag = true;
		else
			WriteRigidContactSurfaces_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteFlexibleContactSurfaces"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteFlexibleContactSurfaces_flag = true;
		else
			WriteFlexibleContactSurfaces_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteForces"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteForces_flag = true;
		else
			WriteForces_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteConstraints"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteConstraints_flag = true;
		else
			WriteConstraints_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteSpecialConstraints"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteSpecialConstraints_flag = true;
		else
			WriteSpecialConstraints_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteContactForces"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteContactForces_flag = true;
		else
			WriteContactForces_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteRenderRigidBodies"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteRenderRigidBodies_flag = true;
		else
			WriteRenderRigidBodies_flag = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteRenderParticles"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteRenderParticles_flag = true;
		else
			WriteRenderParticles_flag = false;
	}
	else
		return false;

	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteSplines"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteSpline_flag = true;
		else
			WriteSpline_flag = false;
	}
	else
	{
		fsetpos(f, &pos);
	}

	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteRenderSplines"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			WriteRenderSpline_flag = true;
		else
			WriteRenderSpline_flag = false;
	}
	else
	{
		fsetpos(f, &pos);
	}
	
	return true;
}

void PostFiles::Write(FILE *f)
{
	fprintf(f, "PostFiles\nMagFactor\t%.6f\nWriteMesh\t%d\nWriteRenderMesh\t%d\nWriteRigidContactSurfaces\t%d\nWriteFlexibleContactSurfaces\t%d\nWriteForces\t%d\nWriteConstraints\t%d\nWriteSpecialConstraints\t%d\nWriteContactForces\t%d\nWriteRenderRigidBodies\t%d\nWriteRenderParticles\t%d\nWriteSplines\t%d\nWriteRenderSplines\t%d\n", mag_factor, WriteMesh_flag, WriteRenderMesh_flag, WriteRigidContactSurfaces_flag, WriteFlexibleContactSurfaces_flag, WriteForces_flag, WriteConstraints_flag, WriteSpecialConstraints_flag, WriteContactForces_flag, WriteRenderRigidBodies_flag, WriteRenderParticles_flag, WriteSpline_flag, WriteRenderSpline_flag);
}

void PostFiles::StartPostFiles(int sol_index)
{
	if (sol_index <= db.number_solutions)
	{
		if (typeid(*db.solution[sol_index - 1]) == typeid(Static) || typeid(*db.solution[sol_index - 1]) == typeid(Dynamic) || typeid(*db.solution[sol_index - 1]) == typeid(ExplicitDynamic))
			StartSinglePartFiles(sol_index);
		if (typeid(*db.solution[sol_index - 1]) == typeid(Modal))
			StartMultiplePartFiles(sol_index);
	}
	else
		StartSinglePartFiles(sol_index);//whole solution
}

void PostFiles::StartSinglePartFiles(int sol_index)
{
	char name[1000];
	//Escrevendo o nome da solução
	char sol_name[100];
	if (sol_index <= db.number_solutions)
		sprintf(sol_name, "solution_%d", sol_index);
	else
		sprintf(sol_name, "whole_solution");
	struct stat info;
	if (WriteMesh_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_mesh.pvd");		//criando arquivo
		f_mesh_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_mesh_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_mesh_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_mesh_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteRenderMesh_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_rendermesh.pvd");
		f_rendermesh_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_rendermesh_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_rendermesh_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_rendermesh_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteFlexibleContactSurfaces_flag == true || WriteRigidContactSurfaces_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_contactsurfaces.pvd");
		f_contactsurfaces_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_contactsurfaces_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_contactsurfaces_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_contactsurfaces_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteSpecialConstraints_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_symbols.pvd");
		f_symbols_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_symbols_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_symbols_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_symbols_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteContactForces_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_contactforces.pvd");
		f_contactforces_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_contactforces_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_contactforces_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_contactforces_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteRenderParticles_flag == true || WriteRenderRigidBodies_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_rb_particles.pvd");
		f_rb_particles_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_rb_particles_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_rb_particles_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_rb_particles_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteConstraints_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_constraints.pvd");
		f_constraints_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_constraints_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_constraints_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_constraints_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteForces_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_forces.pvd");
		f_forces_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_forces_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_forces_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_forces_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteSpline_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_spline.pvd");		//criando arquivo
		f_spline_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_spline_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_spline_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_spline_pvd[sol_index - 1], "\t<Collection>\n");
	}
	if (WriteRenderSpline_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_renderspline.pvd");		//criando arquivo
		f_renderspline_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_renderspline_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_renderspline_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_renderspline_pvd[sol_index - 1], "\t<Collection>\n");
	}
}
void PostFiles::StartMultiplePartFiles(int sol_index)
{
	char name[1000];
	//Escrevendo o nome da solução
	char sol_name[100];
	sprintf(sol_name, "solution_%d", sol_index);

	struct stat info;
	if (WriteMesh_flag == true)
	{
		strcpy(name, db.folder_name);	//pasta do job
		strcat(name, "post/");			//diretório post
		if (stat(name, &info) != 0)		//checa existência do diretório post
			_mkdir(name);				//criando diretório post
		strcat(name, sol_name);			//nome da solution
		strcat(name, "_mesh.pvd");		//criando arquivo
		f_mesh_pvd[sol_index - 1] = fopen(name, "w");
		//Cabeçalho do arquivo VTK
		fprintf(f_mesh_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
		fprintf(f_mesh_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
		fprintf(f_mesh_pvd[sol_index - 1], "\t<Collection>\n");
	}
	//if (WriteRenderMesh_flag == true)
	//{
	//	strcpy(name, data.folder_name);	//pasta do job
	//	strcat(name, "post/");			//diretório post
	//	if (stat(name, &info) != 0)		//checa existência do diretório post
	//		_mkdir(name);				//criando diretório post
	//	strcat(name, sol_name);			//nome da solution
	//	strcat(name, "_rendermesh.pvd");
	//	f_rendermesh_pvd[sol_index - 1] = fopen(name, "w");
	//	//Cabeçalho do arquivo VTK
	//	fprintf(f_rendermesh_pvd[sol_index - 1], "<?xml version=\"1.0\"?>\n");
	//	fprintf(f_rendermesh_pvd[sol_index - 1], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	//	fprintf(f_rendermesh_pvd[sol_index - 1], "\t<Collection>\n");
	//}
}

void PostFiles::UpdateSinglePartPostFiles(int sol_index, double time, int index)
{
	if (WriteMesh_flag == true)
		WriteVTK_XMLBase(sol_index, time, index, 0);
	if (WriteRenderMesh_flag == true)
		WriteVTK_XMLRender(sol_index, time, index, 0);
	if (WriteFlexibleContactSurfaces_flag == true || WriteRigidContactSurfaces_flag == true)
		WriteVTK_XMLContact(sol_index, time, index);
	if (WriteSpecialConstraints_flag == true)
		WriteVTK_XMLSymbols(sol_index, time, index);
	if (WriteContactForces_flag == true)
		WriteVTK_XMLContactForces(sol_index, time, index);
	if (WriteRenderParticles_flag == true || WriteRenderRigidBodies_flag == true)
		WriteVTK_XMLRBParticles(sol_index, time, index);
	if (WriteForces_flag == true)
		WriteVTK_XMLForces(sol_index, time, index);
	if (WriteConstraints_flag == true)
		WriteVTK_XMLConstraints(sol_index, time, index);
	if (WriteSpline_flag == true)
		WriteVTK_XMLSpline(sol_index, time, index);
	if (WriteRenderSpline_flag == true)
		WriteVTK_XMLRenderSpline(sol_index, time, index);
}

void PostFiles::UpdateMultiplePartPostFiles(int sol_index, double time, int index, int part)
{
	if (WriteMesh_flag == true)
		WriteVTK_XMLBase(sol_index, time, index, part);
	/*if (WriteRenderMesh_flag == true)
		WriteVTK_XMLRender(sol_index, time, index, part);*/
}

void PostFiles::EndPostFiles(int sol_index)
{
	if (sol_index <= db.number_solutions)
	{
		if (typeid(*db.solution[sol_index - 1]) == typeid(Static) || typeid(*db.solution[sol_index - 1]) == typeid(Dynamic) || typeid(*db.solution[sol_index - 1]) == typeid(ExplicitDynamic))
			EndSinglePartPostFiles(sol_index);
		if (typeid(*db.solution[sol_index - 1]) == typeid(Modal))
			EndMultiplePartPostFiles(sol_index);
	}
	else
		EndSinglePartPostFiles(sol_index);//whole solution

	
}

void PostFiles::EndSinglePartPostFiles(int sol_index)
{
	if (WriteMesh_flag == true)
	{
		fprintf(f_mesh_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_mesh_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_mesh_pvd[sol_index - 1]);
	}
	if (WriteRenderMesh_flag == true)
	{
		fprintf(f_rendermesh_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_rendermesh_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_rendermesh_pvd[sol_index - 1]);
	}
	if (WriteFlexibleContactSurfaces_flag == true || WriteRigidContactSurfaces_flag == true)
	{
		fprintf(f_contactsurfaces_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_contactsurfaces_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_contactsurfaces_pvd[sol_index - 1]);
	}
	if (WriteSpecialConstraints_flag == true)
	{
		fprintf(f_symbols_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_symbols_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_symbols_pvd[sol_index - 1]);
	}
	if (WriteContactForces_flag == true)
	{
		fprintf(f_contactforces_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_contactforces_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_contactforces_pvd[sol_index - 1]);
	}
	if (WriteRenderParticles_flag == true || WriteRenderRigidBodies_flag == true)
	{
		fprintf(f_rb_particles_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_rb_particles_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_rb_particles_pvd[sol_index - 1]);
	}
	if (WriteForces_flag == true)
	{
		fprintf(f_forces_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_forces_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_forces_pvd[sol_index - 1]);
	}
	if (WriteConstraints_flag == true)
	{
		fprintf(f_constraints_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_constraints_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_constraints_pvd[sol_index - 1]);
	}
	if (WriteSpline_flag == true)
	{
		fprintf(f_spline_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_spline_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_spline_pvd[sol_index - 1]);
	}
	if (WriteRenderSpline_flag == true)
	{
		fprintf(f_renderspline_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_renderspline_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_renderspline_pvd[sol_index - 1]);
	}
}

void PostFiles::EndMultiplePartPostFiles(int sol_index)
{
	if (WriteMesh_flag == true)
	{
		fprintf(f_mesh_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_mesh_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_mesh_pvd[sol_index - 1]);
	}
	/*if (WriteRenderMesh_flag == true)
	{
		fprintf(f_rendermesh_pvd[sol_index - 1], "\t</Collection>\n");
		fprintf(f_rendermesh_pvd[sol_index - 1], "</VTKFile>\n");
		fclose(f_rendermesh_pvd[sol_index - 1]);
	}*/
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLRender(int sol_index, double time, int index, int part)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_rendermesh_");
	strcat(name, number);
	if (part != 0)
	{
		strcat(name, "_part_");
		_itoa(part, number2, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, number2);
	}
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s",analysis,analysis);
	strcat(name, "_rendermesh_");
	strcat(name, number);
	if (part != 0)
	{
		strcat(name, "_part_");
		strcat(name, number2);
	}
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_rendermesh_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"%d\" name=\"Part %d\" file=\"%s\"/>\n", time, part, part, name);
	if (part == 0)
		fprintf(f_rendermesh_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	
	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int ele = 0; ele < db.number_elements; ele++)
	{
		if (typeid(*db.elements[ele]) != typeid(RigidBody_1))
			db.elements[ele]->WriteVTK_XMLRender(f_VTK);
	}
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLRBParticles(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_rb_particles_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_rb_particles_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_rb_particles_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	fprintf(f_rb_particles_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);

	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	if (db.post_files->WriteRenderRigidBodies_flag == true)
	{
		for (int ele = 0; ele < db.number_elements; ele++)
		{
			if (typeid(*db.elements[ele]) == typeid(RigidBody_1))
				db.elements[ele]->WriteVTK_XMLRender(f_VTK);
		}
	}
	if (db.post_files->WriteRenderParticles_flag == true)
	{
		//Conteúdo das Pieces que irão estar no arquivo .vtu
		for (int ele = 0; ele < db.number_particles; ele++)
			db.particles[ele]->WriteVTK_XMLRender(f_VTK);
	}
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}
//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLContact(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_contactsurfaces_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_contactsurfaces_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_contactsurfaces_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	fprintf(f_contactsurfaces_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);

	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int conta = 0; conta < db.number_contacts; conta++)
		db.contacts[conta]->WriteVTK_XMLRender(f_VTK);
	if (WriteRigidContactSurfaces_flag)
	{
		//Conteúdo das Pieces que irão estar no arquivo .vtu
		for (int surf = 0; surf < db.number_boundaries; surf++)
			db.boundaries[surf]->WriteVTK_XMLRender(f_VTK);
	}
	if (WriteFlexibleContactSurfaces_flag)
	{
		//Conteúdo das Geometries que irão estar no arquivo .vtu
		for (int geom = 0; geom < db.number_geometries; geom++)
			db.geometries[geom]->WriteVTK_XMLRender(f_VTK);
	}
	
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLContactForces(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_contactforces_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_contactforces_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_contactforces_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	fprintf(f_contactforces_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	
	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int conta = 0; conta < db.number_contacts; conta++)
		if (db.contacts[conta]->bool_table.GetAt(db.current_solution_number - 1))
			db.contacts[conta]->WriteVTK_XMLForces(f_VTK);
	//General contact search
	if (db.gcs_exist == true)
		db.gcs->WriteVTK_XMLForces(f_VTK);
	
	WriteSingleVertexPart(f_VTK);
	
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLBase(int sol_index, double time, int index, int part)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);
	char number[20];
	char number2[20];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	strcat(name, analysis);
	strcat(name, "_mesh_");
	strcat(name, number);
	if (part != 0)
	{
		strcat(name, "_part_");
		_itoa(part, number2, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, number2);
	}
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");
	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_mesh_");
	strcat(name, number);
	if (part != 0)
	{
		strcat(name, "_part_");
		strcat(name, number2);
	}
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_mesh_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"%d\" name=\"Part %d\" file=\"%s\"/>\n", time, part, part, name);
	//Somente atualiza o pvd do whole solution se part for 0 (ou seja, se não for análise multi-part (modal, por exemplo))
	if (part == 0)
		fprintf(f_mesh_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	
	double tempx, tempy, tempz;
	
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;
	
	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	if (db.number_elements != 0)
	{
		//Opens Piece
		fprintf(f_VTK, "\t\t<Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", db.number_nodes, db.number_elements);
		//Opens Points
		fprintf(f_VTK, "\t\t\t<Points>\n");
		float_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format=\"binary\">\n");
		for (int n = 0; n < db.number_nodes; n++)
		{
			tempx = db.nodes[n]->ref_coordinates[0] + db.post_files->mag_factor*(db.nodes[n]->copy_coordinates[0] - db.nodes[n]->ref_coordinates[0]);
			tempy = db.nodes[n]->ref_coordinates[1] + db.post_files->mag_factor*(db.nodes[n]->copy_coordinates[1] - db.nodes[n]->ref_coordinates[1]);
			tempz = db.nodes[n]->ref_coordinates[2] + db.post_files->mag_factor*(db.nodes[n]->copy_coordinates[2] - db.nodes[n]->ref_coordinates[2]);
			float_vector.push_back((float)(tempx));
			float_vector.push_back((float)(tempy));
			float_vector.push_back((float)(tempz));
		}
		fprintf(f_VTK, encodeData<float>(float_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		//Closes Points
		fprintf(f_VTK, "\t\t\t</Points>\n");

		//Opens Cells
		fprintf(f_VTK, "\t\t\t<Cells>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Int32\" Name = \"connectivity\" format=\"binary\">\n");
		for (int e = 0; e < db.number_elements; e++)
		{
			for (int count = 0; count < db.elements[e]->n_nodes; count++)
				int_vector.push_back(db.elements[e]->nodes[db.elements[e]->VTK_nodes[count]] - 1);
		}
		fprintf(f_VTK, encodeData(int_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Int32\" Name = \"types\" format=\"binary\">\n");
		for (int e = 0; e < db.number_elements; e++)
		{
			int_vector.push_back(db.elements[e]->VTK_type);
			//fprintf(f_VTK, "%d\t", data.elements[e]->VTK_type);
		}
		fprintf(f_VTK, encodeData(int_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Int32\" Name = \"offsets\" format=\"binary\">\n");
		int cur_offset = 0;
		for (int e = 0; e < db.number_elements; e++)
		{
			cur_offset += db.elements[e]->n_nodes;
			//fprintf(f_VTK, "\t\t\t\t\t");
			int_vector.push_back(cur_offset);
			//fprintf(f_VTK, "%d\t", cur_offset);
			//fprintf(f_VTK, "\n");
		}
		fprintf(f_VTK, encodeData(int_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f_VTK, "\t\t\t</Cells>\n");
		if (typeid(*db.solution[db.current_solution_number - 1]) != typeid(Modal))
		{
			//Opens PointData
			fprintf(f_VTK, "\t\t\t<PointData Vectors = \"Displacement Velocity Acceleration\">\n");
			float_vector.clear();
			//Opens DataArray
			fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"Displacement\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
			for (int n = 0; n < db.number_nodes; n++)
			{
				float_vector.push_back((float)(db.nodes[n]->copy_coordinates[0] - db.nodes[n]->ref_coordinates[0]));
				float_vector.push_back((float)(db.nodes[n]->copy_coordinates[1] - db.nodes[n]->ref_coordinates[1]));
				float_vector.push_back((float)(db.nodes[n]->copy_coordinates[2] - db.nodes[n]->ref_coordinates[2]));
			}
			fprintf(f_VTK, encodeData<float>(float_vector).c_str());
			fprintf(f_VTK, "\n");
			//Closes DataArray
			fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
			float_vector.clear();
			//Opens DataArray
			fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"Velocity\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
			for (int n = 0; n < db.number_nodes; n++)
			{
				float_vector.push_back((float)(db.nodes[n]->copy_vel[0]));
				float_vector.push_back((float)(db.nodes[n]->copy_vel[1]));
				float_vector.push_back((float)(db.nodes[n]->copy_vel[2]));
			}
			fprintf(f_VTK, encodeData<float>(float_vector).c_str());
			fprintf(f_VTK, "\n");
			//Closes DataArray
			fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
			float_vector.clear();
			//Opens DataArray
			fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"Acceleration\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
			for (int n = 0; n < db.number_nodes; n++)
			{
				float_vector.push_back((float)(db.nodes[n]->copy_accel[0]));
				float_vector.push_back((float)(db.nodes[n]->copy_accel[1]));
				float_vector.push_back((float)(db.nodes[n]->copy_accel[2]));
			}
			fprintf(f_VTK, encodeData<float>(float_vector).c_str());
			fprintf(f_VTK, "\n");
			//Closes DataArray
			fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
			//Closes PointData
			fprintf(f_VTK, "\t\t\t</PointData>\n");

			//Opens CellData
			fprintf(f_VTK, "\t\t\t<CellData FieldData = \"ElementResults\">\n");
			float_vector.clear();
			//Opens DataArray
			fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"ElementResults\" type = \"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", db.n_element_results);
			for (int e = 0; e < db.number_elements; e++)
				db.elements[e]->WriteVTK_XMLBase(&float_vector);
			fprintf(f_VTK, encodeData<float>(float_vector).c_str());
			fprintf(f_VTK, "\n");
			//Closes DataArray
			fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
			//Closes CellData
			fprintf(f_VTK, "\t\t\t</CellData>\n");
		}
		//Closes Piece
		fprintf(f_VTK, "\t\t</Piece>\n");
	}
	if (db.number_particles != 0)
	{
		//Opens Piece
		fprintf(f_VTK, "\t\t<Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", db.number_particles, db.number_particles);
		//Opens Points
		fprintf(f_VTK, "\t\t\t<Points>\n");
		float_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format=\"binary\">\n");
		for (int i = 0; i < db.number_particles; i++)
		{
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_coordinates[0]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_coordinates[1]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_coordinates[2]));
		}
		fprintf(f_VTK, encodeData<float>(float_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		//Closes Points
		fprintf(f_VTK, "\t\t\t</Points>\n");

		//Opens Cells
		fprintf(f_VTK, "\t\t\t<Cells>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Int32\" Name = \"connectivity\" format=\"binary\">\n");
		for (int e = 0; e < db.number_particles; e++)
		{
			int_vector.push_back(e);
			/*fprintf(f_VTK, "\t\t\t\t\t");
			fprintf(f_VTK, "%d\t", e);
			fprintf(f_VTK, "\n");*/
		}
		fprintf(f_VTK, encodeData(int_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Int32\" Name = \"types\" format=\"binary\">\n");
		for (int e = 0; e < db.number_particles; e++)
		{
			int_vector.push_back(1);
			//fprintf(f_VTK, "\t\t\t\t\t");
			//fprintf(f_VTK, "%d\t", 1);//VTK vertex
			//fprintf(f_VTK, "\n");
		}
		fprintf(f_VTK, encodeData(int_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray type = \"Int32\" Name = \"offsets\" format=\"binary\">\n");
		int cur_offset = 0;
		for (int e = 0; e < db.number_particles; e++)
		{
			cur_offset++;
			int_vector.push_back(cur_offset);
		}
		fprintf(f_VTK, encodeData(int_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f_VTK, "\t\t\t</Cells>\n");

		//Opens PointData
		fprintf(f_VTK, "\t\t\t<PointData Vectors = \"Displacement Velocity Acceleration\" Scalars = \"Diameter\">\n");
		float_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"Displacement\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int i = 0; i < db.number_particles; i++)
		{
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_coordinates[0] - db.nodes[db.particles[i]->node - 1]->ref_coordinates[0]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_coordinates[1] - db.nodes[db.particles[i]->node - 1]->ref_coordinates[1]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_coordinates[2] - db.nodes[db.particles[i]->node - 1]->ref_coordinates[2]));
		}
		fprintf(f_VTK, encodeData<float>(float_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		float_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"Velocity\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int i = 0; i < db.number_particles; i++)
		{
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_vel[0]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_vel[1]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_vel[2]));
		}
		fprintf(f_VTK, encodeData<float>(float_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");
		float_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"Acceleration\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int i = 0; i < db.number_particles; i++)
		{
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_accel[0]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_accel[1]));
			float_vector.push_back((float)(db.nodes[db.particles[i]->node - 1]->copy_accel[2]));
		}
		fprintf(f_VTK, encodeData<float>(float_vector).c_str());
		fprintf(f_VTK, "\n");
		float_vector.clear();
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");

		//Diameter
		float_vector.clear();
		//Opens DataArray
		fprintf(f_VTK, "\t\t\t\t<DataArray Name = \"Diameter\" type = \"Float32\" NumberOfComponents= \"1\" format=\"binary\">\n");
		for (int i = 0; i < db.number_particles; i++)
		{
			if (typeid(*db.particles[i]) == typeid(Sphere))
			{
				//Ponteiro para o Sphere
				Sphere* ptr = static_cast<Sphere*>(db.particles[i]);
				float_vector.push_back((float)(ptr->radius * 2));
			}
			else
			{
				float_vector.push_back((float)(0.0));
			}
		}
		fprintf(f_VTK, encodeData<float>(float_vector).c_str());
		fprintf(f_VTK, "\n");
		//Closes DataArray
		fprintf(f_VTK, "\t\t\t\t</DataArray>\n");

		//Closes PointData
		fprintf(f_VTK, "\t\t\t</PointData>\n");
		//Closes Piece
		fprintf(f_VTK, "\t\t</Piece>\n");
	}
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLSpline(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_spline_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_spline_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd	
	fprintf(f_spline_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"%d\" name=\"Part %d\" file=\"%s\"/>\n", time, 0, 0, name);
	//Somente atualiza o pvd do whole solution se part for 0 (ou seja, se não for análise multi-part (modal, por exemplo))
	fprintf(f_spline_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);

	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int ele = 0; ele < db.number_splines; ele++)
	{
		db.splines[ele]->WriteVTK_XML_SplineMesh(f_VTK);
	}
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLRenderSpline(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_renderspline_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_renderspline_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd	
	fprintf(f_renderspline_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"%d\" name=\"Part %d\" file=\"%s\"/>\n", time, 0, 0, name);
	//Somente atualiza o pvd do whole solution se part for 0 (ou seja, se não for análise multi-part (modal, por exemplo))
	fprintf(f_renderspline_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);

	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int ele = 0; ele < db.number_splines; ele++)
	{
		db.splines[ele]->WriteVTK_XML_SplineRender(f_VTK);
	}
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLSymbols(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_symbols_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_symbols_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_symbols_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	fprintf(f_symbols_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	
	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int i = 0; i < db.number_special_constraints; i++)
		db.special_constraints[i]->WriteVTK_XMLRender(f_VTK);
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLConstraints(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_constraints_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_constraints_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_constraints_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	fprintf(f_constraints_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);

	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int i = 0; i < db.number_constraints; i++)
		db.constraints[i]->WriteVTK_XML(f_VTK);
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

//Gerando arquivo para o Paraview (VTK) - formato XML
void PostFiles::WriteVTK_XMLForces(int sol_index, double time, int index)
{
	FILE *f_VTK;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_forces_");
	strcat(name, number);
	strcat(name, ".vtu");
	f_VTK = fopen(name, "w");

	//Nome do arquivo do pvd
	sprintf(name, "./%s/%s", analysis, analysis);
	strcat(name, "_forces_");
	strcat(name, number);
	strcat(name, ".vtu");
	//Atualiza pvd
	fprintf(f_forces_pvd[sol_index - 1], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	fprintf(f_forces_pvd[db.number_solutions], "\t\t<DataSet timestep=\"%.10f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, name);
	
	//Cabeçalho do arquivo VTK
	fprintf(f_VTK, "<?xml version=\"1.0\"?>\n");
	fprintf(f_VTK, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f_VTK, "<!--INPUT: %s - TIME = %.10f-->\n", db.file_name, time);
	fprintf(f_VTK, "\t<UnstructuredGrid>\n");
	//Conteúdo das Pieces que irão estar no arquivo .vtu
	for (int i = 0; i < db.number_loads; i++)
		db.loads[i]->WriteVTK_XML(f_VTK);
	fprintf(f_VTK, "\t</UnstructuredGrid>\n");
	fprintf(f_VTK, "</VTKFile>\n");
	fclose(f_VTK);
}

void PostFiles::WriteConfigurationResults(int sol_index, double time, int index)
{
	FILE *f;
	char name[1000];
	strcpy(name, db.folder_name);	//pasta do job
	struct stat info;
	strcat(name, "post/");			//diretório post
	if (stat(name, &info) != 0)		//checa existência do diretório post
		_mkdir(name);				//criando diretório post
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", sol_index);
	strcat(name, analysis);
	strcat(name, "/");
	_mkdir(name);	//cria a pasta com o tipo de análise em questão
	char number[20];
	char number2[50];
	_itoa(index - 1, number, 10);	//Converte o número inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
	sprintf(number2, "%lf", time);	//Converte o time para char
	strcat(name, analysis);
	strcat(name, "_configuration_");
	strcat(name, number);
	strcat(name, ".txt");
	f = fopen(name, "w");

	fprintf(f, "RESULTS - NODES\n\n");
	for (int i = 0; i < db.number_nodes; i++)
		db.nodes[i]->WriteResults(f);

	fprintf(f, "\n\nRESULTS - ELEMENTS\n\n");
	for (int i = 0; i < db.number_elements; i++)
		db.elements[i]->WriteResults(f);

	fclose(f);
}

void PostFiles::WriteSingleVertexPart(FILE *f)
{
	//Creates a single null contact force at origin - to avoid bugs in paraview visualization
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	float_vector.clear();
	float_vector.push_back(0.0);
	float_vector.push_back(0.0);
	float_vector.push_back(0.0);
	std::vector<int> int_vector;
	//Opens Piece
	fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 1, 1);
	//Opens Points
	fprintf(f, "\t\t\t<Points>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Points
	fprintf(f, "\t\t\t</Points>\n");
	//Opens Cells
	fprintf(f, "\t\t\t<Cells>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	int_vector.clear();
	int_vector.push_back(0);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	int_vector.clear();
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	int_vector.push_back(1);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Cells
	fprintf(f, "\t\t\t</Cells>\n");
	//Opens PointData
	fprintf(f, "\t\t\t<PointData Vectors = \"ContactForces\">\n");
	fprintf(f, "\t\t\t\t<DataArray Name = \"Normal\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	fprintf(f, "\t\t\t\t<DataArray Name = \"Friction\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes PointData
	fprintf(f, "\t\t\t</PointData>\n");
	//Closes Piece
	fprintf(f, "\t\t</Piece>\n");
}