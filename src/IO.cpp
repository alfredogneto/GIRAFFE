#include "IO.h"
#include <string.h>
#include <sys/stat.h>
#include "PlatformCompat.h"
#define PI 3.1415926535897932384626433832795

#include "Node.h"
#include "SuperNode.h"
#include "PipeSection.h"
#include "RigidBodyData.h"
#include "CoordinateSystem.h"
#include "Environment.h"

#include "Monitor.h"
#include "PostFiles.h"
#include "SolverOptions.h"

#include "NodeSet.h"
#include "SurfaceSet.h"
#include "ElementSet.h"
#include "SuperNodeSet.h"

#include "InitialCondition.h"
#include "ConvergenceCriteria.h"

#include "AnalyticalSurface.h"
#include "Plane.h"
#include "LineRegion.h"
#include "SurfaceRegion.h"

//Elements
#include "Element.h"
#include "Beam_1.h"
#include "Pipe_1.h"
#include "Shell_1.h"
#include "Solid_1.h"
#include "SpringDashpot_1.h"
#include "Mass_1.h"
#include "RigidBody_1.h"
#include "Truss_1.h"
#include "TwoNodeConnector_1.h"

//Materials
#include "Material.h"
#include "Hooke.h"
#include "ElasticPlasticIsoHardening.h"
#include "Orthotropic.h"

//Point & Curves
#include "Point.h"
#include "ArcCirc.h"

//Surfaces
#include "Surface.h"
#include "RigidTriangularSurface_1.h"
#include "RigidOscillatorySurface_1.h"
#include "FlexibleTriangularSurface_2.h"
#include "FlexibleSECylinder_1.h"
#include "FlexibleArcExtrusion_1.h"
#include "RigidArcRevolution_1.h"
#include "RigidNURBS_1.h"

//Geometries
#include "Geometry.h"
#include "SECylinder.h"
#include "ArcExtrusion.h"
#include "ArcRevolution.h"

//BodyGeometry
#include "BodyGeometry.h"

//Splines
#include "Spline.h"

//Contacts
#include "Contact.h"
#include "LRLR.h"
#include "GeneralPLR.h"
#include "NSSS.h"
#include "SSSS.h"
#include "SPSP.h"

//GeneralContactSearch
#include "GeneralContactSearch.h"

//Sections transversais
#include "Section.h"
#include "SecGeneral.h"
#include "SecRectangle.h"
#include "SecSuperEllipse.h"
#include "SecTube.h"
#include "SecUserDefined.h"
#include "SecHelicalFiber.h"

//Shell sections
#include "ShellSection.h"
#include "ShellSectionHomogeneous.h"
#include "ShellSectionComposite.h"

//Special constraints
#include "SpecialConstraint.h"
#include "SameDisplacement.h"
#include "Hinge.h"
#include "UniversalJoint.h"
#include "SameRotation.h"
#include "RigidNodeSet.h"
#include "TranslationalJoint.h"
#include "NodalConstraintDOF.h"

//Particles
#include "Particle.h"
#include "Sphere.h"
#include "Polyhedron.h"
#include "NURBSParticle.h"
#include "VEMPolyhedron.h"

//Boundaries
#include "Boundary.h"
#include "STLBoundary.h"

//Surface pairs (for general contact search)
#include "SurfacePairGeneralContact.h"
#include "RigidTriangularFace_RigidTriangularFace.h"
#include "FlexibleTriangularFace_FlexibleTriangularFace.h"
#include "FlexibleTriangularFace_RigidTriangularFace.h"

//Section details
#include "SectionDetails.h"
#include "SolidSection.h"
#include "MultiCellSection.h"

//AerodynamicData
#include "AerodynamicData.h"
#include "BEM.h"

//Solutions
#include "Solution.h"
#include "Static.h"
#include "Dynamic.h"
#include "Modal.h"
#include "ConcomitantSolution.h"
#include "ExplicitDynamic.h"

//Loads
#include "Load.h"
#include "NodalLoad.h"
#include "NodalFollowerLoad.h"
#include "ShellLoad.h"
#include "PipeLoad.h"
#include "SuperNodalLoad.h"

//Displacements prescritos
#include "Displacement.h"
#include "NodalDisplacement.h"
#include "DisplacementField.h"

//Constraints
#include "Constraint.h"
#include "NodalConstraint.h"
#include "SuperNodalConstraint.h"

//PSY
#include "PSYCoupling.h"

//CADData
#include "CADData.h"
#include "STLSurface.h"
#include "NURBSSurface.h"

//ContactInterfaces
#include "ContactInterface.h"
#include "Interface_0.h"
#include "Interface_1.h"

//BoundingVolumes
#include "BoundingVolume.h"
#include "BoundingSphere.h"
#include "BoundingCylinder.h"
#include "BoundingTriangularBox.h"
#include "BoundingBoxAxesAligned.h"
#include "BoundingBoxAxesOriented.h"

#include "Encoding.h"
#include "ExecutionData.h"
#include "ConfigurationSave.h"

IO::IO(void)
{
}

IO::~IO(void)
{
}
bool IO::ReadFile(int argc, char* argv[])
{
	//Reading the input file name
	char name_input [1000];
	bool readOK = false;
	FILE *f = NULL;
	printf(" ____________________________________________________________ \n");
	printf("|                                                            |\n");
	printf("|                          GIRAFFE                           |\n");
	printf("|                                                            |\n");
	printf("|  Generic Interface Readily Accessible for Finite Elements  |\n");
	printf("|                                                            |\n");
	printf("|               University of Sao Paulo - Brazil             |\n");
	printf("|                                                            |\n");
	printf("|____________________________________________________________|\n\n"); 
	while (readOK == false)
	{
		if (argc <= 1)
		{
			printf("Enter the name of the input file:\n");
			scanf("%s", name_input);
		}
		else
		{
			strcpy(name_input, argv[1]);
		}
		strcpy(db.file_name, name_input);
		strcat(db.file_name, ".inp");
		strcpy(db.folder_name, "./");
		strcat(db.folder_name, name_input);
		strcat(db.folder_name, "/");
		char name[1000];
		strcpy(name, db.folder_name);
		strcat(name, db.file_name);
		printf("\n");
		//tries to read the same location of the executable file of Giraffe
		f = fopen(name, "r");
		//if no folder/file is found, the second try is on the Giraffe instalation directory, inside the folder inputs
		if (f == NULL)
		{
			//saves directory of the current loged user on windows
			char* dir_install;
			//saves the name of the required environment variable
			char variable[20];
			strcpy(variable, "GIRAFFE_INSTALL");
			//reads the environment variable
			dir_install = std::getenv(variable);
			strcpy(db.folder_name, dir_install);
			strcat(db.folder_name, "/inputs/");
			strcat(db.folder_name, name_input);
			strcat(db.folder_name, "/");
			strcpy(name, db.folder_name);
			strcat(name, db.file_name);
			f = fopen(name, "r");
			//last try is on the public directory, folder /Documents/Giraffe/
			if (f == NULL)
			{
				//saves directory of the current loged user (PUBLIC on Windows, HOME on macOS/Linux)
				char* dir_user;
				//reads the environment variable
				dir_user = std::getenv("PUBLIC");
				if (dir_user == NULL)
					dir_user = std::getenv("HOME");
				if (dir_user == NULL)
				{
					printf("Error reading the input file.\n");
					return false;
				}
				strcpy(db.folder_name, dir_user);
				strcat(db.folder_name, "/Documents/Giraffe/");
				//checks if there is a folder called "Giraffe" in this directory
				struct stat info;
				if (stat(db.folder_name, &info) != 0)
					GIRAFFE_MKDIR(db.folder_name);
				strcat(db.folder_name, name_input);
				strcat(db.folder_name, "/");
				strcpy(name, db.folder_name);
				strcat(name, db.file_name);
				f = fopen(name, "r");
				if (f == NULL)
				{
					printf("Error reading the input file.\n");
					return false;
				}
				else
					readOK = true;
			}
			else
				readOK = true;
		}
		else
			readOK = true;
	}
	char s[1000];
	char last_key[1000];
	bool read = true;
	bool EOF_key_read = false;
	while (fscanf(f, "%s", s) != EOF && EOF_key_read == false)
	{
		read = true;
		///////////////////////////Reading of the Solution////////////////////////////
		if (!strcmp(s, "SolutionSteps") && read == true)
		{
			if (!ReadSolutions(f))
				return false;//Error
			read = false;
			db.solution_exist = true;
		}
		///////////////////////////Reading of the nodes////////////////////////////////
		if (!strcmp(s, "Nodes") && read == true)
		{
			if (!ReadNodes(f))	
				return false;//Error
			read = false;
			db.nodes_exist = true;
		}
		///////////////////////////Reading of the super nodes/////////////////////////
		if (!strcmp(s, "SuperNodes") && read == true)
		{
			if (!ReadSuperNodes(f))
				return false;//Error
			read = false;
			db.super_nodes_exist = true;
		}
		///////////////////////////Reading of the points//////////////////////////////
		if (!strcmp(s, "Points") && read == true)
		{
			if (!ReadPoints(f))
				return false;//Error
			read = false;
			db.points_exist = true;
		}
		///////////////////////////Reading of the arcos////////////////////////////////
		if (!strcmp(s, "Arcs") && read == true)
		{
			if (!ReadArcs(f))
				return false;//Error
			read = false;
			db.arcs_exist = true;
		}
		///////////////////////////Reading of the elements////////////////////////////
		if (!strcmp(s, "Elements") && read == true)
		{
			if (!ReadElements(f))
				return false;//Error
			read = false;
			db.elements_exist = true;
		}
		///////////////////////////Reading of the particles////////////////////////////
		if (!strcmp(s, "Particles") && read == true)
		{
			if (!ReadParticles(f))
				return false;//Error
			read = false;
			db.particles_exist = true;
		}
		///////////////////////////Reading of the InitialConditions////////////////////////
		if (!strcmp(s, "InitialConditions") && read == true)
		{
			if (!ReadInitialConditions(f))
				return false;//Error
			read = false;
			db.IC_exist = true;
		}
		///////////////////////////Reading of the materials////////////////////////////
		if (!strcmp(s, "Materials") && read == true)
		{
			if (!ReadMaterials(f))
				return false;//Error
			read = false;
			db.materials_exist = true;
		}
		///////////////////////////Reading of the sections//////////////////////////////
		if (!strcmp(s, "Sections") && read == true)
		{
			if (!ReadSections(f))
				return false;//Error
			read = false;
			db.sections_exist = true;
		}
		///////////////////////////Reading of the sections of tubos/////////////////////
		if (!strcmp(s, "PipeSections") && read == true)
		{
			if (!ReadPipeSections(f))
				return false;//Error
			read = false;
			db.pipe_sections_exist = true;
		}
		///////////////////////////Reading of the sections of shells/////////////////////
		if (!strcmp(s, "ShellSections") && read == true)
		{
			if (!ReadShellSections(f))
				return false;//Error
			read = false;
			db.shell_sections_exist = true;
		}
		///////////////////////////Reading of the CS//////////////////////////////////
		if (!strcmp(s, "CoordinateSystems") && read == true)
		{
			if (!ReadCoordinateSystems(f))
				return false;//Error
			read = false;
			db.CS_exist = true;
		}
		///////////////////////////Reading of the propriedades of rigid body/////////////////////
		if (!strcmp(s, "RigidBodyData") && read == true)
		{
			if (!ReadRigidBodyData(f))
				return false;//Error
			read = false;
			db.RB_data_exist = true;
		}
		///////////////////////////Reading Environment/////////////////////////////
		if (!strcmp(s, "Environment") && read == true)
		{
			if (!ReadEnvironment(f))
				return false;//Error
			read = false;
			db.environment_exist = true;
		}
		///////////////////////////Reading of the Monitor////////////////////////////////
		if (!strcmp(s, "Monitor") && read == true)
		{
			if (!ReadMonitor(f))
				return false;//Error
			read = false;
			db.monitor_exist = true;
		}
		///////////////////////////Reading of the SolverOptions////////////////////////////
		if (!strcmp(s, "SolverOptions") && read == true)
		{
			if (!ReadSolverOptions(f))
				return false;//Error
			read = false;
			db.solver_options_exist = true;
		}
		///////////////////////////Reading of the Analytical Surfaces//////////////////
		if (!strcmp(s, "AnalyticalSurfaces"))
		{
			if (!ReadAnalyticalSurfaces(f))
				return false;//Error
			read = false;
			db.analytical_surfaces_exist = true;
		}
		///////////////////////////Reading of the Surfaces//////////////////
		if (!strcmp(s, "Surfaces"))
		{
			if (!ReadSurfaces(f))
				return false;//Error
			read = false;
			db.surfaces_exist = true;
		}
		///////////////////////////Reading of the Splines////////////////////////////////
		if (!strcmp(s, "Splines") && read == true)
		{
			if (!ReadSplines(f))
				return false;//Error
			read = false;
			db.splines_exist = true;
		}
		///////////////////////////Reading of the LineRegions//////////////////////////
		if (!strcmp(s, "LineRegions"))
		{
			if (!ReadLineRegions(f))
				return false;//Error
			read = false;
			db.line_regions_exist = true;
		}
		///////////////////////////Reading of the SurfaceRegions//////////////////////////
		if (!strcmp(s, "SurfaceRegions"))
		{
			if (!ReadSurfaceRegions(f))
				return false;//Error
			read = false;
			db.surface_regions_exist = true;
		}
		///////////////////////////Reading of the Contacts/////////////////////////////
		if (!strcmp(s, "Contacts"))
		{
			if (!ReadContacts(f))
				return false;//Error
			read = false;
			db.contacts_exist = true;
		}
		///////////////////////////Reading of the NodeSets////////////////////////////
		if (!strcmp(s, "NodeSets") && read == true)
		{
			if (!ReadNodeSets(f))
				return false;//Error
			read = false;
			db.node_sets_exist = true;
		}
		///////////////////////////Reading of the SuperNodeSets////////////////////////////
		if (!strcmp(s, "SuperNodeSets") && read == true)
		{
			if (!ReadSuperNodeSets(f))
				return false;//Error
			read = false;
			db.super_node_sets_exist = true;
		}
		///////////////////////////Reading of the SurfaceSets////////////////////////////
		if (!strcmp(s, "SurfaceSets") && read == true)
		{
			if (!ReadSurfaceSets(f))
				return false;//Error
			read = false;
			db.surface_sets_exist = true;
		}
		///////////////////////////Reading of the ElementSets////////////////////////////
		if (!strcmp(s, "ElementSets") && read == true)
		{
			if (!ReadElementSets(f))
				return false;//Error
			read = false;
			db.element_sets_exist = true;
		}
		///////////////////////////Reading of the Loads////////////////////////////////
		if (!strcmp(s, "Loads") && read == true)
		{
			if (!ReadLoads(f))
				return false;//Error
			read = false;
			db.loads_exist = true;
		}
		///////////////////////////Reading of the Displacements////////////////////////
		if (!strcmp(s, "Displacements") && read == true)
		{
			if (!ReadDisplacements(f))
				return false;//Error
			read = false;
			db.displacements_exist = true;
		}
		/////////////////////////Reading of the Constraints/////////////////////////////
		if (!strcmp(s, "Constraints"))
		{
			if (!ReadConstraints(f))
				return false;//Error
			read = false;
			db.constraints_exist = true;
		}
		/////////////////////////Reading of the Special Constraints///////////////////////
		if (!strcmp(s, "SpecialConstraints"))
		{
			if (!ReadSpecialConstraints(f))
				return false;//Error
			read = false;
			db.special_constraints_exist = true;
		}
		/////////////////////////Reading of the Section Details///////////////////////
		if (!strcmp(s, "SectionDetails"))
		{
			if (!ReadSectionDetails(f))
				return false;//Error
			read = false;
			db.section_details_exist = true;
		}
		///////////////////////////Data reading aerodinamicos//////////////////
		if (!strcmp(s, "AerodynamicData") && read == true)
		{
			if (!ReadAerodynamicData(f))
				return false;//Error
			read = false;
			db.aerodynamic_data_exist = true;
		}
		///////////////////////////Reading BEM/////////////////////////////
		if (!strcmp(s, "BEM") && read == true)
		{
			if (!ReadBEMData(f))
				return false;//Error
			read = false;
			db.bem_exist = true;
		}
		///////////////////////////Reading GeneralContactSearch/////////////////////////////
		if (!strcmp(s, "GeneralContactSearch") && read == true)
		{
			if (!ReadGeneralContactSearch(f))
				return false;//Error
			read = false;
			db.gcs_exist = true;
		}
		///////////////////////////Reading ParticlePack/////////////////////////////
		if (!strcmp(s, "ConfigurationSave") && read == true)
		{
			if (!ReadConfigurationSave(f))
				return false;//Error
			read = false;
			db.config_save_exist = true;
		}
		///////////////////////////Reading ConcomitantModalSolution///////////////
		if (!strcmp(s, "ConcomitantSolution") && read == true)
		{
			if (!ReadConcomitantSolution(f))
				return false;//Error
			read = false;
			db.concomitant_solution_exist = true;
		} 
		///////////////////////////Reading of the ConvergenceCriteria////////////////////
		if (!strcmp(s, "ConvergenceCriteria"))
		{
			if (!ReadConvergenceCriteria(f))
				return false;//Error
			read = false;
		}
		///////////////////////////Reading of the PostFiles///////////////////////////////
		if (!strcmp(s, "PostFiles") && read == true)
		{
			if (!ReadPostFiles(f))
				return false;//Error
			read = false;
		}
		///////////////////////////Reading of the execution data///////////////////////
		if (!strcmp(s, "ExecutionData") && read == true)
		{
			if (!ReadExecutionData(f))
				return false;//Error
			read = false;
		}
		///////////////////////////Reading of the PSYCoupling////////////////////////////////
		if (!strcmp(s, "PSYCoupling") && read == true)
		{
			if (!ReadPSYCoupling(f))
				return false;//Error
			read = false;
			db.psy_coupling_exist = true;
		}

		///////////////////////////Reading of the CADs////////////////////////////
		if (!strcmp(s, "CADData") && read == true)
		{
			if (!ReadCADData(f))
				return false;//Error
			read = false;
			db.cad_data_exist = true;
		}
		///////////////////////////Reading of the contact interfaces////////////////////////////
		if (!strcmp(s, "ContactInterfaces") && read == true)
		{
			if (!ReadContactInterfaces(f))
				return false;//Error
			read = false;
			db.contactinterfaces_exist = true;
		}

		///////////////////////////Reading of the boundaries////////////////////////////
		if (!strcmp(s, "Boundaries") && read == true)
		{
			if (!ReadBoundaries(f))
				return false;//Error
			read = false;
			db.boundaries_exist = true;
		}
		///////////////////////////Reading of the geometries////////////////////////////
		if (!strcmp(s, "Geometries") && read == true)
		{
			if (!ReadGeometries(f))
				return false;//Error
			read = false;
			db.geometries_exist = true;
		}
		///////////////////////////Reading of the BodyGeometries////////////////////////////
		if (!strcmp(s, "BodyGeometries") && read == true)
		{
			if (!ReadBodyGeometries(f))
				return false;//Error
			read = false;
			db.body_geometries_exist = true;
		}
		///////////////////////////End of the reading - keyword that indicates end////
		if (!strcmp(s, "EOF"))
		{
			EOF_key_read = true;
			read = false;
		}
		
		//If read == true in this point, trata-if of keyword not conhecida
		if (read == true)
		{
			//Calls reading of possible comment in the input file - returns true if there was reading of comment - false if not there was
			read = ReadComment(f, s, 10000);
			if (read == false)
			{
				printf("Error reading a keyword after %s block.\n", last_key);
				return false;//Error
			}
			
		}
		else
			strcpy(last_key, s);//Copies the keyword for saber qual was the last lida
		
	}
	fclose(f);
	//Abertura of the streaming for save output in file of text
	strcpy(db.name, db.folder_name);
	strcat(db.name, "simulation_report.txt");
	db.console_output = fopen(db.name, "w");
	return true;
}

//Reads comments - returns the stream in the point after reading of comment
bool ReadComment(FILE *f, char* s, int dim_char)
{
	bool read = false;
	///////////////////////////Comment///////////////////////////////////////
	if (s[0] == '/' && s[1] == '/')//linha of comment
	{
		//Reading until find the caracter of end of linha
		while (s[0] != '\n')
			fscanf(f, "%c", s);
		read = true;
	}
	if (s[0] == '/' && s[1] == '*')//bloco of comment
	{
		//Reading until find the end of the comment
		char c1[1];
		char c2[1];
		bool exit = false;
		//Search through the end of the bloco in two fases
		//1) In the próprio text already read saved in "s"
		for (int i = 2; i < 999; i++)
		{
			if (s[i] == '*' && s[i + 1] == '/')
				exit = true;
		}
		//2) Outside of the próprio text already read. Search of caracter in caracter until find the end of the comment
		while (exit == false)
		{
			//Reading of the first caracter
			fscanf(f, "%c", c1);
			//If encontrou the first caracter of the output of the comment - tries check the second caracter for confirmar output
			if (c1[0] == '*')
			{
				//Saves the position (stream) - before of the reading of the second caractere
				fpos_t pos;
				fgetpos(f, &pos);
				fscanf(f, "%c", c2);
				if (c2[0] == '/')
					exit = true;
				else
					fsetpos(f, &pos);
			}
		}
		read = true;
	}
	return read;
}

//Tries read comments. Returns the stream in the point in that the next reading not and a comment
void TryComment(FILE *f)
{
	bool comment = true;//flag that indicates that leu comment in the last attempt
	char s[10000];
	while (comment == true)
	{
		//Saves the position (stream) - before of the reading
		fpos_t pos;
		fgetpos(f, &pos);
		//Reading
		fscanf(f, "%s", s);
		comment = ReadComment(f, s, 10000);
		if (comment == false)
			fsetpos(f, &pos);//turn the position original of the stream - not leu comment
	}
	
}

void WriteDOFTable(FILE *f)
{
	fprintf(f, "Node\tUX\tUY\tUZ\tROTX\tROTY\tROTZ\n");
	for (int node = 0; node < db.number_nodes; node++)
	{
		fprintf(f, "%d\t", node + 1);
		for (int k = 0; k < db.number_GLs_node; k++)
		{
			if (db.nodes[node]->GLs[k] > 0)
				fprintf(f, "%d\t", db.nodes[node]->GLs[k]);
			else
				fprintf(f, "fixed\t");
		}
		fprintf(f, "\n");
	}
}

void IO::WriteFile()
{
	//Fechamento of the streaming for save output in file of text
	fclose(db.console_output);
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "output.inp");
	FILE *f = fopen(name, "w");
	if (db.solution_exist == true)
		WriteSolutions(f);
	if (db.nodes_exist == true)
		WriteNodes(f);
	if (db.super_nodes_exist == true)
		WriteSuperNodes(f);
	if (db.points_exist == true)
		WritePoints(f);
	if (db.arcs_exist == true)
		WriteArcs(f);
	if (db.elements_exist == true)
		WriteElements(f);
	if (db.particles_exist == true)
		WriteParticles(f);
	if (db.IC_exist == true)
		WriteInitialConditions(f);
	if (db.materials_exist == true)
		WriteMaterials(f);
	if (db.sections_exist == true)
		WriteSections(f);
	if (db.pipe_sections_exist == true)
		WritePipeSections(f);
	if (db.shell_sections_exist == true)
		WriteShellSections(f);
	if (db.CS_exist == true)
		WriteCoordinateSystems(f);
	if (db.RB_data_exist == true)
		WriteRigidBodyData(f);
	if (db.environment_exist == true)
		WriteEnvironment(f);
	if (db.monitor_exist == true)
		WriteMonitor(f);
	if (db.solver_options_exist == true)
		WriteSolverOptions(f);
	if (db.analytical_surfaces_exist == true)
		WriteAnalyticalSurfaces(f);
	if (db.surfaces_exist == true)
		WriteSurfaces(f);
	if (db.splines_exist == true)
		WriteSplines(f);
	if (db.line_regions_exist == true)
		WriteLineRegions(f);
	if (db.surface_regions_exist == true)
		WriteSurfaceRegions(f);
	if (db.contacts_exist == true)
		WriteContacts(f);
	if (db.node_sets_exist == true)
		WriteNodeSets(f);
	if (db.super_node_sets_exist == true)
		WriteSuperNodeSets(f);
	if (db.surface_sets_exist == true)
		WriteSurfaceSets(f);
	if (db.element_sets_exist == true)
		WriteElementSets(f);
	if (db.loads_exist == true)
		WriteLoads(f);
	if (db.displacements_exist == true)
		WriteDisplacements(f);
	if (db.constraints_exist == true)
		WriteConstraints(f);
	if (db.special_constraints_exist == true)
		WriteSpecialConstraints(f);
	if (db.section_details_exist == true)
		WriteSectionDetails(f);
	if (db.aerodynamic_data_exist == true)
		WriteAerodynamicData(f);
	if (db.cad_data_exist == true)
		WriteCADData(f);
	if (db.contactinterfaces_exist == true)
		WriteContactInterfaces(f);
	if (db.geometries_exist == true)
		WriteGeometries(f);
	if (db.body_geometries_exist == true)
		WriteBodyGeometries(f);
	if (db.bem_exist == true)
		WriteBEMData(f);
	if (db.gcs_exist == true)
		WriteGeneralContactSearch(f);
	if (db.concomitant_solution_exist == true)
		WriteConcomitantSolution(f);
	WriteConvergenceCriteria(f);
	WritePostFiles(f);
	WriteExecutionData(f);
	fclose(f);
}


bool IO::ReadSolutions(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	db.number_solutions = atoi(s);
	db.solution = new Solution*[db.number_solutions];	//Allocation of the vector
	bool solution_OK = false;
	for (int i = 0; i < db.number_solutions; i++)
	{
		TryComment(f);
		solution_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the solution "Static"
		if (!strcmp(s, "Static"))
		{
			db.solution[i] = new Static();
			solution_OK = true;
		}
		//Allocation of the solution "Dynamic"
		if (!strcmp(s, "Dynamic"))
		{
			db.solution[i] = new Dynamic();
			solution_OK = true;
		}
		//Allocation of the solution "ExplicitDynamic"
		if (!strcmp(s, "ExplicitDynamic"))
		{
			db.solution[i] = new ExplicitDynamic();
			solution_OK = true;
		}
		//Allocation of the solution "Modal"
		if (!strcmp(s, "Modal"))
		{
			db.solution[i] = new Modal();
			solution_OK = true;
		}
		if (solution_OK == true)
		{
			if (!db.solution[i]->Read(f))
			{
				printf("Error reading Solution %d.\n", i + 1);
				db.number_solutions = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Solution %d or %d.\n", i, i + 1);
			db.number_solutions = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	db.EvaluateStartEndTimes();	//Fills data of start and end time of the solutions - and importante do this filling here and not in the PreCalc(), since the function of check errors and called in the sequência imediata, before of the PreCalc()

	return true;
}

bool IO::ReadNodes(FILE *f)
{
	char s[1000];
	int n_nodes = 0;
	fscanf(f, "%s", s);
	n_nodes = atoi(s);
	db.number_nodes = n_nodes;		//sets in the database
	db.nodes = new Node*[n_nodes];	//Allocation of the vector
	for (int i = 0; i < n_nodes; i++)
	{
		db.nodes[i] = new Node(db.number_GLs_node);		//Allocation of each node
		TryComment(f);
		if (!db.nodes[i]->Read(f))						//Data reading of the node
		{
			printf("Error reading Node %d.\n",i+1);
			db.number_nodes = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadPoints(FILE *f)
{
	char s[1000];
	int n_points = 0;
	fscanf(f, "%s", s);
	n_points = atoi(s);
	db.number_points = n_points;		//sets in the database
	db.points = new Point*[n_points];	//Allocation of the vector
	for (int i = 0; i < n_points; i++)
	{
		db.points[i] = new Point();		//Allocation of each point
		TryComment(f);
		if (!db.points[i]->Read(f))		//Data reading of the point
		{
			printf("Error reading Point %d.\n", i + 1);
			db.number_points = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadArcs(FILE *f)
{
	char s[1000];
	int n_arcs = 0;
	fscanf(f, "%s", s);
	n_arcs = atoi(s);
	db.number_arcs = n_arcs;		//sets in the database
	db.arcs = new ArcCirc*[n_arcs];	//Allocation of the vector
	for (int i = 0; i < n_arcs; i++)
	{
		db.arcs[i] = new ArcCirc();		//Allocation of each point
		TryComment(f);
		if (!db.arcs[i]->Read(f))		//Data reading of the point
		{
			printf("Error reading Arc %d.\n", i + 1);
			db.number_arcs = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadElements(FILE *f)
{
	char s[1000];
	int n_elements = 0;
	fscanf(f, "%s", s);
	n_elements = atoi(s);
	db.number_elements = n_elements;			//sets in the database
	db.elements = new Element*[n_elements];	//Allocation of the vector
	bool element_OK = false;
	for (int i = 0; i < n_elements; i++)
	{
		TryComment(f);
		element_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the element "Beam_1"
		if (!strcmp(s, "Beam_1"))
		{
			db.elements[i] = new Beam_1();
			element_OK = true;
		}
		//Allocation of the element "Pipe_1"
		if (!strcmp(s, "Pipe_1"))
		{
			db.elements[i] = new Pipe_1();
			element_OK = true;
		}
		//Allocation of the element "Shell_1"
		if (!strcmp(s, "Shell_1"))
		{
			db.elements[i] = new Shell_1();
			element_OK = true;
		}
		//Allocation of the element "Solid_1"
		if (!strcmp(s, "Solid_1"))
		{
			db.elements[i] = new Solid_1();
			element_OK = true;
		}
		//Allocation of the element "SpringDashpot_1"
		if (!strcmp(s, "SpringDashpot_1"))
		{
			db.elements[i] = new SpringDashpot_1();
			element_OK = true;
		}
		//Allocation of the element "Mass_1"
		if (!strcmp(s, "Mass_1"))
		{
			db.elements[i] = new Mass_1();
			element_OK = true;
		}
		//Allocation of the element "RigidBody_1"
		if (!strcmp(s, "RigidBody_1"))
		{
			db.elements[i] = new RigidBody_1();
			element_OK = true;
		}
		//Allocation of the element "Truss_1"
		if (!strcmp(s, "Truss_1"))
		{
			db.elements[i] = new Truss_1();
			element_OK = true;
		}
		//Allocation of the element "TwoNodeConnector_1"
		if (!strcmp(s, "TwoNodeConnector_1"))
		{
			db.elements[i] = new TwoNodeConnector_1();
			element_OK = true;
		}
		if (element_OK == true)
		{
			if (!db.elements[i]->Read(f))
			{
				printf("Error reading Element %d.\n", i + 1);
				db.number_elements = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Element %d or %d.\n", i, i + 1);
			db.number_elements = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadParticles(FILE *f)
{
	char s[1000];
	int n_particles = 0;
	fscanf(f, "%s", s);
	n_particles = atoi(s);
	db.number_particles = n_particles;		//sets in the database
	db.particles = new Particle*[n_particles];	//Allocation of the vector
	bool particle_OK = false;
	for (int i = 0; i < n_particles; i++)
	{
		TryComment(f);
		particle_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the particle "Sphere"
		if (!strcmp(s, "Sphere"))
		{
			db.particles[i] = new Sphere();
			particle_OK = true;
		}

		//Allocation of the particle "Polyhedron"
		if (!strcmp(s, "Polyhedron"))
		{
			db.particles[i] = new Polyhedron();
			particle_OK = true;
		}

		//Allocation of the particle "NURBS"
		if (!strcmp(s, "NURBSParticle"))
		{
			db.particles[i] = new NURBSParticle();
			particle_OK = true;
		}

		//Allocation of the particle "VEMPolyhedron"
		if (!strcmp(s, "VEMPolyhedron"))
		{
			db.particles[i] = new VEMPolyhedron();
			particle_OK = true;
		}
		
		if (particle_OK == true)
		{
			
			if (!db.particles[i]->Read(f))
			{
				printf("Error reading Particle %d.\n", i + 1);
				db.number_particles = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Particle %d or %d.\n", i, i + 1);
			db.number_particles = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadInitialConditions(FILE *f)
{
	char s[1000];
	int n_IC = 0;
	fscanf(f, "%s", s);
	n_IC = atoi(s);
	db.number_IC = n_IC;	//sets in the database

	db.IC = new InitialCondition*[n_IC];	//Allocation of the vector
	for (int i = 0; i < n_IC; i++)
	{
		db.IC[i] = new InitialCondition();		//Allocation of each initial condition
		TryComment(f);
		if (!db.IC[i]->Read(f))
		{
			printf("Error reading Initial Condition %d.\n", i + 1);
			db.number_IC = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadMaterials(FILE *f)
{
	char s[1000];
	int n_materials = 0;
	fscanf(f, "%s", s);
	n_materials = atoi(s);
	db.number_materials = n_materials;			//sets in the database
	db.materials = new Material*[n_materials];	//Allocation of the vector
	bool materialOK = false;
	for (int i = 0; i < n_materials; i++)
	{
		TryComment(f);
		materialOK = false;
		fscanf(f, "%s", s);
		//Allocation of the material "Hooke"
		if (!strcmp(s, "Hooke"))
		{
			db.materials[i] = new Hooke();
			materialOK = true;
		}
		//Allocation of the material "ElasticPlasticIsoHardening"
		if (!strcmp(s, "ElasticPlasticIsoHardening"))
		{
			db.materials[i] = new ElasticPlasticIsoHardening();
			materialOK = true;
		}
		//Allocation of the material "Orthotropic"
		if (!strcmp(s, "Orthotropic"))
		{
			db.materials[i] = new Orthotropic();
			materialOK = true;
		}
		
		if (materialOK == true)
		{
			
			if (!db.materials[i]->Read(f))
			{
				printf("Error reading Material %d.\n", i + 1);
				db.number_materials = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Material %d or %d.\n",i, i + 1);
			db.number_materials = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}
bool IO::ReadSections(FILE *f)
{
	char s[1000];
	int n_sec = 0;
	fscanf(f, "%s", s);
	n_sec = atoi(s);
	db.number_sections = n_sec;	//sets in the database
	bool sectionOK = false;
	db.sections = new Section*[n_sec];	//Allocation of the vector
	for (int i = 0; i < n_sec; i++)
	{
		TryComment(f);
		sectionOK = false;
		fscanf(f, "%s", s);
		//Allocation of the section "General"
		if (!strcmp(s, "General"))
		{
			db.sections[i] = new SecGeneral();		//Allocation of the section
			sectionOK = true;
		}
		//Allocation of the section "Rectangle"
		if (!strcmp(s, "Rectangle"))
		{
			db.sections[i] = new SecRectangle();		//Allocation of the section
			sectionOK = true;
		}	
		//Allocation of the section "SuperEllipse"
		if (!strcmp(s, "SuperEllipse"))
		{
			db.sections[i] = new SecSuperEllipse();	//Allocation of the section
			sectionOK = true;
		}
		//Allocation of the section "SuperEllipse"
		if (!strcmp(s, "Tube"))
		{
			db.sections[i] = new SecTube();			//Allocation of the section
			sectionOK = true;
		}
		//Allocation of the section "UserDefined"
		if (!strcmp(s, "UserDefined"))
		{
			db.sections[i] = new SecUserDefined();		//Allocation of the section
			sectionOK = true;
		}
		//Allocation of the section "HelicalFiber"
		if (!strcmp(s, "HelicalFiber"))
		{
			db.sections[i] = new SecHelicalFiber();		//Allocation of the section
			sectionOK = true;
		}
		if (sectionOK == true)
		{
			
			if (!db.sections[i]->Read(f))
			{
				printf("Error reading Section %d.\n", i + 1);
				db.number_sections = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Section %d or %d.\n",i, i + 1);
			db.number_sections = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}
bool IO::ReadPipeSections(FILE *f)
{
	char s[1000];
	int n_pipe_sec = 0;
	fscanf(f, "%s", s);
	n_pipe_sec = atoi(s);
	db.number_pipe_sections = n_pipe_sec;	//sets in the database

	db.pipe_sections = new PipeSection*[n_pipe_sec];	//Allocation of the vector
	for (int i = 0; i < n_pipe_sec; i++)
	{
		db.pipe_sections[i] = new PipeSection();		//Allocation of each Pipe Section
		TryComment(f);
		if (!db.pipe_sections[i]->Read(f))
		{
			printf("Error reading Pipe Section %d.\n", i + 1);
			db.number_pipe_sections = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}
bool IO::ReadShellSections(FILE *f)
{
	char s[1000];
	int n_shell_sec = 0;
	fscanf(f, "%s", s);
	n_shell_sec = atoi(s);
	db.number_shell_sections = n_shell_sec;	//sets in the database
	bool section_OK = false;
	db.shell_sections = new ShellSection*[n_shell_sec];	//Allocation of the vector
	for (int i = 0; i < n_shell_sec; i++)
	{
		TryComment(f);
		section_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the element "Homogeneous"
		if (!strcmp(s, "Homogeneous"))
		{
			db.shell_sections[i] = new ShellSectionHomogeneous();
			section_OK = true;
		}
		//Allocation of the element "Composite"
		if (!strcmp(s, "Composite"))
		{
			db.shell_sections[i] = new ShellSectionComposite();
			section_OK = true;
		}
		
		if (section_OK == true)
		{
			if (!db.shell_sections[i]->Read(f))
			{
				printf("Error reading Shell Section %d.\n", i + 1);
				db.number_shell_sections = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Shell Section %d or %d.\n", i, i + 1);
			db.number_shell_sections = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}
bool IO::ReadCoordinateSystems(FILE *f)
{
	char s[1000];
	int n_CS = 0;
	fscanf(f, "%s", s);
	n_CS = atoi(s);
	db.number_CS = n_CS;	//sets in the database
	db.CS = new CoordinateSystem*[n_CS];	//Allocation of the vector
	for (int i = 0; i < n_CS; i++)
	{
		db.CS[i] = new CoordinateSystem();		//Allocation of each CS
		//Data reading of the CS
		TryComment(f);
		if (!db.CS[i]->Read(f))
		{
			printf("Error reading CS %d.\n", i + 1);
			db.number_CS = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadRigidBodyData(FILE *f)
{
	char s[1000];
	int n_RBD = 0;
	fscanf(f, "%s", s);
	n_RBD = atoi(s);
	db.number_RB_data = n_RBD;	//sets in the database
	db.RB_data = new RigidBodyData*[n_RBD];	//Allocation of the vector
	for (int i = 0; i < n_RBD; i++)
	{
		db.RB_data[i] = new RigidBodyData();		//Allocation of each CS
		TryComment(f);
		//Data reading of the RBD
		if (!db.RB_data[i]->Read(f))
		{
			printf("Error reading Rigid Body Data %d.\n", i + 1);
			db.number_RB_data = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadEnvironment(FILE *f)
{
	db.environment = new Environment();
	TryComment(f);
	if (!db.environment->Read(f))
	{
		printf("Error reading Environment.\n");
		return false;
	}
	return true;
}

bool IO::ReadMonitor(FILE *f)
{
	db.monitor = new Monitor();
	TryComment(f);
	if (!db.monitor->Read(f))
	{
		printf("Error reading Monitor.\n");
		return false;
	}
	return true;
}

bool IO::ReadSolverOptions(FILE *f)
{
	TryComment(f);
	if (!db.solver_options->Read(f))
	{
		printf("Error reading Solver Options.\n");
		return false;
	}
	return true;
}

bool IO::ReadAnalyticalSurfaces(FILE *f)
{
	char s[1000];
	int n_analytical = 0;
	fscanf(f, "%s", s);
	n_analytical = atoi(s);
	db.number_analytical_surfaces = n_analytical;	//sets in the database
	bool analyticalOK = false;
	db.analytical_surfaces = new AnalyticalSurface*[n_analytical];	//Allocation of the vector
	for (int i = 0; i < n_analytical; i++)
	{
		TryComment(f);
		analyticalOK = false;
		fscanf(f, "%s", s);
		//Allocation of the "Plane"
		if (!strcmp(s, "Plane"))
		{
			db.analytical_surfaces[i] = new Plane();		//Allocation of each object
			analyticalOK = true;
		}
			
		//Data reading of the object
		if (analyticalOK == true)
		{
			if (!db.analytical_surfaces[i]->Read(f))
			{
				printf("Error reading Analytical Surface %d.\n", i + 1);
				db.number_analytical_surfaces = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Analytical Surface %d or %d.\n",i, i + 1);
			db.number_analytical_surfaces = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
		
	}
	return true;
}

bool IO::ReadSurfaces(FILE *f)
{
	char s[1000];
	int n_surf = 0;
	fscanf(f, "%s", s);
	n_surf = atoi(s);
	db.number_surfaces = n_surf;	//sets in the database
	bool surfOK = false;
	db.surfaces = new Surface*[n_surf];	//Allocation of the vector
	for (int i = 0; i < n_surf; i++)
	{
		TryComment(f);
		surfOK = false;
		fscanf(f, "%s", s);
		//Allocation of the "RigidTriangularSurface_1"
		if (!strcmp(s, "RigidTriangularSurface_1"))
		{
			db.surfaces[i] = new RigidTriangularSurface_1();		//Allocation of each object
			surfOK = true;
		}
		//Allocation of the "RigidOscillatorySurface_1"
		if (!strcmp(s, "RigidOscillatorySurface_1"))
		{
			db.surfaces[i] = new RigidOscillatorySurface_1();		//Allocation of each object
			surfOK = true;
		}
		//Allocation of the "FlexibleTriangularSurface_2"
		if (!strcmp(s, "FlexibleTriangularSurface_2"))
		{
			db.surfaces[i] = new FlexibleTriangularSurface_2();		//Allocation of each object
			surfOK = true;
		}
		//Allocation of the "FlexibleSECylinder_1"
		if (!strcmp(s, "FlexibleSECylinder_1"))
		{
			db.surfaces[i] = new FlexibleSECylinder_1();		//Allocation of each object
			surfOK = true;
		}
		//Allocation of the "FlexibleArcExtrusion_1"
		if (!strcmp(s, "FlexibleArcExtrusion_1"))
		{
			db.surfaces[i] = new FlexibleArcExtrusion_1();		//Allocation of each object
			surfOK = true;
		}

		//Allocation of the "RigidArcRevolution_1"
		if (!strcmp(s, "RigidArcRevolution_1"))
		{
			db.surfaces[i] = new RigidArcRevolution_1();		//Allocation of each object
			surfOK = true;
		}

		//Allocation of the "RigidNURBS_1"
		if (!strcmp(s, "RigidNURBS_1"))
		{
			db.surfaces[i] = new RigidNURBS_1();		//Allocation of each object
			surfOK = true;
		}

		//Data reading of the object
		if (surfOK == true)
		{
			if (!db.surfaces[i]->Read(f))
			{
				printf("Error reading Surface %d.\n", i + 1);
				db.number_surfaces = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Surface %d or %d.\n", i, i + 1);
			db.number_surfaces = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}

	}
	return true;
}

bool IO::ReadSplines(FILE *f)
{
	char s[1000];
	int number_splines = 0;
	fscanf(f, "%s", s);
	number_splines = atoi(s);
	db.number_splines = number_splines;		//sets in the database
	db.splines = new Spline*[number_splines];	//Allocation of the vector
	for (int i = 0; i < number_splines; i++)
	{
		db.splines[i] = new Spline();		//Allocation of each node
		TryComment(f);
		if (!db.splines[i]->Read(f))						//Data reading of the node
		{
			printf("Error reading Spline %d.\n", i + 1);
			db.number_splines = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}


bool IO::ReadLineRegions(FILE *f)
{
	char s[1000];
	int n_line_r = 0;
	fscanf(f, "%s", s);
	n_line_r = atoi(s);
	db.number_line_regions = n_line_r;	//sets in the database

	db.line_regions = new LineRegion*[n_line_r];	//Allocation of the vector
	for (int i = 0; i < n_line_r; i++)
	{
		db.line_regions[i] = new LineRegion();		//Allocation of each object
		TryComment(f);
		//Data reading of the object
		if (!db.line_regions[i]->Read(f))
		{
			printf("Error reading Line Region %d.\n", i + 1);
			db.number_line_regions = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadSurfaceRegions(FILE *f)
{
	char s[1000];
	int n_surf_r = 0;
	fscanf(f, "%s", s);
	n_surf_r = atoi(s);
	db.number_surface_regions = n_surf_r;	//sets in the database

	db.surface_regions = new SurfaceRegion*[n_surf_r];	//Allocation of the vector
	for (int i = 0; i < n_surf_r; i++)
	{
		db.surface_regions[i] = new SurfaceRegion();		//Allocation of each object
		TryComment(f);
		//Data reading of the object
		if (!db.surface_regions[i]->Read(f))
		{
			printf("Error reading Surface Region %d.\n", i + 1);
			db.number_surface_regions = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadContacts(FILE *f)
{
	char s[1000];
	int n_contacts = 0;
	fscanf(f, "%s", s);
	n_contacts = atoi(s);
	db.number_contacts = n_contacts;	//sets in the database
	bool contactOK = false;
	db.contacts = new Contact*[n_contacts];	//Allocation of the vector
	for (int i = 0; i < n_contacts; i++)
	{
		TryComment(f);
		contactOK = false;
		fscanf(f, "%s", s);
		//Allocation of the contact "LRLR"
		if (!strcmp(s, "LRLR"))
		{
			db.contacts[i] = new LRLR();		//Allocation of each object
			contactOK = true;
		}
		//Allocation of the contact "GeneralPLR"
		if (!strcmp(s, "GeneralPLR"))
		{
			db.contacts[i] = new GeneralPLR();	//Allocation of each object
			contactOK = true;
		}
		//Allocation of the contact "NSSS"
		if (!strcmp(s, "NSSS"))
		{
			db.contacts[i] = new NSSS();	//Allocation of each object
			contactOK = true;
		}
		//Allocation of the contact "SSSS"
		if (!strcmp(s, "SSSS"))
		{
			db.contacts[i] = new SSSS();	//Allocation of each object
			contactOK = true;
		}
		//Allocation of the contact "SPSP"
		if (!strcmp(s, "SPSP"))
		{
			db.contacts[i] = new SPSP();	//Allocation of each object
			contactOK = true;
		}
	
		if (contactOK == true)
		{
			//Data reading of the object
			if (!db.contacts[i]->Read(f))
			{
				printf("Error reading Contact %d.\n", i + 1);
				db.number_contacts = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Contact %d or %d.\n",i, i + 1);
			db.number_contacts = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
		
	}
	return true;
}


bool IO::ReadNodeSets(FILE *f)
{
	char s[1000];
	int n_n_sets = 0;
	fscanf(f, "%s", s);
	n_n_sets = atoi(s);
	db.number_node_sets = n_n_sets;	//sets in the database

	db.node_sets = new NodeSet*[n_n_sets];	//Allocation of the vector
	for (int i = 0; i < n_n_sets; i++)
	{
		db.node_sets[i] = new NodeSet();		//Allocation of each object
		TryComment(f);
		//Data reading of the object
		if (!db.node_sets[i]->Read(f))
		{
			printf("Error reading Node Set %d.\n", i + 1);
			db.number_node_sets = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadSuperNodeSets(FILE *f)
{
	char s[1000];
	int n_sn_sets = 0;
	fscanf(f, "%s", s);
	n_sn_sets = atoi(s);
	db.number_super_node_sets = n_sn_sets;	//sets in the database

	db.super_node_sets = new SuperNodeSet*[n_sn_sets];	//Allocation of the vector
	for (int i = 0; i < n_sn_sets; i++)
	{
		db.super_node_sets[i] = new SuperNodeSet();		//Allocation of each object
		TryComment(f);
		//Data reading of the object
		if (!db.super_node_sets[i]->Read(f))
		{
			printf("Error reading Super Node Set %d.\n", i + 1);
			db.number_node_sets = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadSurfaceSets(FILE *f)
{
	char s[1000];
	int n_s_sets = 0;
	fscanf(f, "%s", s);
	n_s_sets = atoi(s);
	db.number_surface_sets = n_s_sets;	//sets in the database
	db.surface_sets = new SurfaceSet*[n_s_sets];	//Allocation of the vector
	for (int i = 0; i < n_s_sets; i++)
	{
		db.surface_sets[i] = new SurfaceSet();		//Allocation of each object
		TryComment(f);
		//Data reading of the object
		if (!db.surface_sets[i]->Read(f))
		{
			printf("Error reading Surface Set %d.\n", i + 1);
			db.number_surface_sets = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadElementSets(FILE *f)
{
	char s[1000];
	int n_e_sets = 0;
	fscanf(f, "%s", s);
	n_e_sets = atoi(s);
	db.number_element_sets = n_e_sets;	//sets in the database
	db.element_sets = new ElementSet*[n_e_sets];	//Allocation of the vector
	for (int i = 0; i < n_e_sets; i++)
	{
		db.element_sets[i] = new ElementSet();		//Allocation of each object
		TryComment(f);
		//Data reading of the object
		if (!db.element_sets[i]->Read(f))
		{
			printf("Error reading Element Set %d.\n", i + 1);
			db.number_element_sets = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadLoads(FILE *f)
{
	char s[1000];
	int n_loads = 0;
	fscanf(f, "%s", s);
	n_loads = atoi(s);
	db.number_loads = n_loads;									//sets in the database
	db.loads = new Load*[n_loads];								//Allocation of the vector
	bool load_OK = false;
	for (int i = 0; i < n_loads; i++)
	{
		TryComment(f);
		load_OK = false;
		fscanf(f, "%s", s);

		//Allocation of the constraint of the type "NodalLoad"
		if (!strcmp(s, "NodalLoad"))
		{
			db.loads[i] = new NodalLoad();
			load_OK = true;
		}
		//Allocation of the constraint of the type "NodalFollowerLoad"
		if (!strcmp(s, "NodalFollowerLoad"))
		{
			db.loads[i] = new NodalFollowerLoad();
			load_OK = true;
		}
		//Allocation of the constraint of the type "PipeLoad"
		if (!strcmp(s, "PipeLoad"))
		{
			db.loads[i] = new PipeLoad();
			load_OK = true;
		}
		//Allocation of the constraint of the type "ShellLoad"
		if (!strcmp(s, "ShellLoad"))
		{
			db.loads[i] = new ShellLoad();
			load_OK = true;
		}
		//Allocation of the constraint of the type "SuperNodalLoad"
		if (!strcmp(s, "SuperNodalLoad"))
		{
			db.loads[i] = new SuperNodalLoad();
			load_OK = true;
		}

		if (load_OK == true)
		{
			if (!db.loads[i]->Read(f))
			{
				printf("Error reading Load %d.\n", i + 1);
				db.number_loads = i + 1;			//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Load %d or %d.\n", i, i + 1);
			db.number_loads = i;					//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadDisplacements(FILE *f)
{
	char s[1000];
	int n_disp = 0;
	fscanf(f, "%s", s);
	n_disp = atoi(s);
	db.number_displacements = n_disp;								//sets in the database
	db.displacements = new Displacement*[n_disp];					//Allocation of the vector
	bool disp_OK = false;
	for (int i = 0; i < n_disp; i++)
	{
		TryComment(f);
		disp_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the constraint of the type "NodalDisplacement"
		if (!strcmp(s, "NodalDisplacement"))
		{
			db.displacements[i] = new NodalDisplacement();
			disp_OK = true;
		}

		//Allocation of the constraint of the type "NodalDisplacement"
		if (!strcmp(s, "DisplacementField"))
		{
			db.displacements[i] = new DisplacementField();
			disp_OK = true;
		}
		
		if (disp_OK == true)
		{
			if (!db.displacements[i]->Read(f))
			{
				printf("Error reading Displacement %d.\n", i + 1);
				db.number_displacements = i + 1;			//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Displacement %d or %d.\n", i, i + 1);
			db.number_displacements = i;					//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}


bool IO::ReadConstraints(FILE *f)
{
	char s[1000];
	int n_const = 0;
	fscanf(f, "%s", s);
	n_const = atoi(s);
	db.number_constraints = n_const;					//sets in the database
	db.constraints = new Constraint*[n_const];		//Allocation of the vector
	bool c_OK = false;
	for (int i = 0; i < n_const; i++)
	{
		TryComment(f);
		c_OK = false;
		fscanf(f, "%s", s);

		//Allocation of the constraint of the type "NodalConstraint"
		if (!strcmp(s, "NodalConstraint"))
		{
			db.constraints[i] = new NodalConstraint();
			c_OK = true;
		}
		//Allocation of the constraint of the type "SuperNodalConstraint"
		if (!strcmp(s, "SuperNodalConstraint"))
		{
			db.constraints[i] = new SuperNodalConstraint();
			c_OK = true;
		}

		if (c_OK == true)
		{
			if (!db.constraints[i]->Read(f))
			{
				printf("Error reading Constraint %d.\n", i + 1);
				db.number_constraints = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Constraint %d or %d.\n", i, i + 1);
			db.number_constraints = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadSpecialConstraints(FILE *f)
{
	char s[1000];
	int n_special_c = 0;
	fscanf(f, "%s", s);
	n_special_c = atoi(s);
	db.number_special_constraints = n_special_c;					//sets in the database
	db.special_constraints = new SpecialConstraint*[n_special_c];	//Allocation of the vector
	bool special_c_OK = false;
	for (int i = 0; i < n_special_c; i++)
	{
		TryComment(f);
		special_c_OK = false;
		fscanf(f, "%s", s);

		//Allocation of the constraint of the type "SameDisplacement"
		if (!strcmp(s, "SameDisplacement"))
		{
			db.special_constraints[i] = new SameDisplacement();
			special_c_OK = true;
		}
		//Allocation of the constraint of the type "Hinge"
		if (!strcmp(s, "HingeJoint"))
		{
			db.special_constraints[i] = new Hinge();
			special_c_OK = true;
		}
		//Allocation of the constraint of the type "UniversalJoint"
		if (!strcmp(s, "UniversalJoint"))
		{
			db.special_constraints[i] = new UniversalJoint();
			special_c_OK = true;
		}
		//Allocation of the constraint of the type "SameRotation"
		if (!strcmp(s, "SameRotation"))
		{
			db.special_constraints[i] = new SameRotation();
			special_c_OK = true;
		}
		//Allocation of the constraint of the type "SameRotation"
		if (!strcmp(s, "RigidNodeSet"))
		{
			db.special_constraints[i] = new RigidNodeSet();
			special_c_OK = true;
		}
		//Allocation of the constraint of the type "TranslationalJoint"
		if (!strcmp(s, "TranslationalJoint"))
		{
			db.special_constraints[i] = new TranslationalJoint();
			special_c_OK = true;
		}
		//Allocation of the constraint of the type "NodalConstraintDOF"
		if (!strcmp(s, "NodalConstraintDOF"))
		{
			db.special_constraints[i] = new NodalConstraintDOF();
			special_c_OK = true;
		}

		if (special_c_OK == true)
		{

			if (!db.special_constraints[i]->Read(f))
			{
				printf("Error reading Special Constraint %d.\n", i + 1);
				db.number_special_constraints = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Special Constraint %d or %d.\n", i, i + 1);
			db.number_special_constraints = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadSectionDetails(FILE *f)
{
	char s[1000];
	int n_sec_details = 0;
	fscanf(f, "%s", s);
	n_sec_details = atoi(s);
	db.number_section_details = n_sec_details;	//sets in the database
	bool read_OK = false;
	db.section_details = new SectionDetails*[n_sec_details];	//Allocation of the vector
	for (int i = 0; i < n_sec_details; i++)
	{
		TryComment(f);
		fscanf(f, "%s", s);
		if (!strcmp(s, "SolidSection"))
		{
			db.section_details[i] = new SolidSection();			//Allocation of each SectionDetails
			read_OK = true;
		}
		if (!strcmp(s, "MultiCellSection"))
		{
			db.section_details[i] = new MultiCellSection();		//Allocation of each SectionDetails
			read_OK = true;
		}
		if (read_OK == true)
		{

			//Data reading of the SectionDetails				
			if (!db.section_details[i]->Read(f))
			{
				printf("Error reading Section Details %d.\n", i + 1);
				db.number_section_details = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Section Details %d or %d.\n", i, i + 1);
			db.number_section_details = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}

	}
	return true;

}

bool IO::ReadAerodynamicData(FILE *f)
{
	char s[1000];
	int n_aero = 0;
	fscanf(f, "%s", s);
	n_aero = atoi(s);
	db.number_aerodynamicdata = n_aero;					//sets in the database
	db.aerodynamic_data = new AerodynamicData*[n_aero];	//Allocation of the vector
	for (int i = 0; i < n_aero; i++)
	{
		db.aerodynamic_data[i] = new AerodynamicData();		//Allocation of each point
		TryComment(f);
		if (!db.aerodynamic_data[i]->Read(f))		//Data reading of the aerodynamicdata
		{
			printf("Error reading AerodynamicData %d.\n", i + 1);
			db.number_aerodynamicdata = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadCADData(FILE *f)
{
	char s[1000];
	int n_CAD = 0;
	fscanf(f, "%s", s);
	n_CAD = atoi(s);
	db.number_cad_data = n_CAD;				//sets in the database
	db.cad_data = new CADData*[n_CAD];		//Allocation of the vector
	bool CAD_OK = false;
	for (int i = 0; i < n_CAD; i++)
	{
		TryComment(f);
		CAD_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the element "STLSurface"
		if (!strcmp(s, "STLSurface"))
		{
			db.cad_data[i] = new STLSurface();
			CAD_OK = true;
		}
		//Allocation of the element "NURBSSurface"
		if (!strcmp(s, "NURBSSurface"))
		{
			db.cad_data[i] = new NURBSSurface();
			CAD_OK = true;
		}
		
		if (CAD_OK == true)
		{
			if (!db.cad_data[i]->Read(f))
			{
				printf("Error reading CADData %d.\n", i + 1);
				db.number_cad_data = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading CADData %d or %d.\n", i, i + 1);
			db.number_cad_data = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadContactInterfaces(FILE *f)
{
	char s[1000];
	int n_interfaces = 0;
	fscanf(f, "%s", s);
	n_interfaces = atoi(s);
	db.number_contactinterfaces = n_interfaces;						//Sets in the database
	db.contactinterfaces = new ContactInterface*[n_interfaces];		//Allocation of the vector
	bool Interface_OK = false;
	for (int i = 0; i < n_interfaces; i++)
	{
		TryComment(f);
		Interface_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the interface "Interface_1"
		if (!strcmp(s, "Interface_1"))
		{
			db.contactinterfaces[i] = new Interface_1();
			Interface_OK = true;
		}
		//Allocation of the interface "Interface_0"
		if (!strcmp(s, "Interface_0"))
		{
			db.contactinterfaces[i] = new Interface_0();
			Interface_OK = true;
		}
		if (Interface_OK == true)
		{
			if (!db.contactinterfaces[i]->Read(f))
			{
				printf("Error reading ContactInterface %d.\n", i + 1);
				db.number_contactinterfaces = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading ContactInterface %d or %d.\n", i, i + 1);
			db.number_contactinterfaces = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadGeometries(FILE *f)
{
	char s[1000];
	int n_geometries = 0;
	fscanf(f, "%s", s);
	n_geometries = atoi(s);
	db.number_geometries = n_geometries;							//Sets in the database
	db.geometries = new Geometry*[n_geometries];					//Allocation of the vector
	bool GEOM_OK = false;
	for (int i = 0; i < n_geometries; i++)
	{
		TryComment(f);
		GEOM_OK = false;
		fscanf(f, "%s", s);
		//Allocation of the interface "SECylinder"
		if (!strcmp(s, "SECylinder"))
		{
			db.geometries[i] = new SECylinder();
			GEOM_OK = true;
		}
		//Allocation of the interface "ArcExtrusion"
		if (!strcmp(s, "ArcExtrusion"))
		{
			db.geometries[i] = new ArcExtrusion();
			GEOM_OK = true;
		}
		//Allocation of the interface "ArcRevolution"
		if (!strcmp(s, "ArcRevolution"))
		{
			db.geometries[i] = new ArcRevolution();
			GEOM_OK = true;
		}

		if (GEOM_OK == true)
		{
			if (!db.geometries[i]->Read(f))
			{
				printf("Error reading Geometry %d.\n", i + 1);
				db.number_geometries = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Geometry %d or %d.\n", i, i + 1);
			db.number_geometries = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadBodyGeometries(FILE *f)
{
	char s[1000];
	int n_body = 0;
	fscanf(f, "%s", s);
	n_body = atoi(s);
	db.number_body_geometries = n_body;					//sets in the database
	db.body_geometries = new BodyGeometry*[n_body];		//Allocation of the vector
	for (int i = 0; i < n_body; i++)
	{
		db.body_geometries[i] = new BodyGeometry();		//Allocation of each body
		TryComment(f);
		if (!db.body_geometries[i]->Read(f))			//Data reading
		{
			printf("Error reading BodyGeometry %d.\n", i + 1);
			db.number_body_geometries = i + 1;					//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

bool IO::ReadBoundaries(FILE *f)
{
	char s[1000];
	int n_bound = 0;
	fscanf(f, "%s", s);
	n_bound = atoi(s);
	db.number_boundaries = n_bound;						//Sets in the database
	db.boundaries = new Boundary*[n_bound];				//Allocation of the vector
	bool boundary_OK = false;
	for (int i = 0; i < n_bound; i++)
	{
		TryComment(f);
		boundary_OK = false;
		fscanf(f, "%s", s);
		//Allocation of "STLBoundary"
		if (!strcmp(s, "STLBoundary"))
		{
			db.boundaries[i] = new STLBoundary();
			boundary_OK = true;
		}
		if (boundary_OK == true)
		{
			if (!db.boundaries[i]->Read(f))
			{
				printf("Error reading Boundary %d.\n", i + 1);
				db.number_boundaries = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
				return false;
			}
		}
		else
		{
			printf("Error reading Boundary %d or %d.\n", i, i + 1);
			db.number_boundaries = i;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}


bool IO::ReadBEMData(FILE *f)
{
	db.bem = new BEM();
	TryComment(f);
	if (!db.bem->Read(f))
	{
		printf("Error reading BEM.\n");
		return false;
	}
	return true;
}

bool IO::ReadGeneralContactSearch(FILE *f)
{
	db.gcs = new GeneralContactSearch();
	TryComment(f);
	if (!db.gcs->Read(f))
	{
		printf("Error reading GeneralContactSearch.\n");
		return false;
	}
	return true;
}

bool IO::ReadConfigurationSave(FILE *f)
{
	db.config_save = new ConfigurationSave();
	TryComment(f);
	if (!db.config_save->Read(f))
	{
		printf("Error reading ConfigurationSave.\n");
		return false;
	}
	return true;
}

bool IO::ReadConcomitantSolution(FILE *f)
{
	db.concomitant_solution = new ConcomitantSolution();
	TryComment(f);
	if (!db.concomitant_solution->Read(f))
	{
		printf("Error reading ConcomitantSolution.\n");
		return false;
	}
	return true;
}

bool IO::ReadConvergenceCriteria(FILE *f)
{
	TryComment(f);
	if (!db.conv_criteria->Read(f))
	{
		printf("Error reading Convergence Criteria.\n");
		return false;
	}
	return true;
}

bool IO::ReadPostFiles(FILE *f)
{
	TryComment(f);
	if (!db.post_files->Read(f))
	{
		printf("Error reading Post Files.\n");
		return false;
	}
	return true;
}

bool IO::ReadPSYCoupling(FILE *f)
{
	db.psy_coupling = new PSYCoupling();
	TryComment(f);
	if (!db.psy_coupling->Read(f))
	{
		printf("Error reading PSY Coupling.\n");
		return false;
	}
	return true;
}

bool IO::ReadExecutionData(FILE* f)
{
	TryComment(f);
	if (!db.execution_data->Read(f))
	{
		printf("Error reading Execution Data.\n");
		return false;
	}
	return true;
}
bool IO::ReadSuperNodes(FILE *f)
{
	char s[1000];
	int n_sn = 0;
	fscanf(f, "%s", s);
	n_sn = atoi(s);
	db.number_super_nodes = n_sn;			//sets in the database
	db.super_nodes = new SuperNode*[n_sn];	//Allocation of the vector
	for (int i = 0; i < n_sn; i++)
	{
		db.super_nodes[i] = new SuperNode();			//Allocation of each super node
		TryComment(f);
		if (!db.super_nodes[i]->Read(f))						//Data reading of the super node
		{
			printf("Error reading Super Node %d.\n", i + 1);
			db.number_super_nodes = i + 1;	//changes the number of allocated instances - avoids an error in the destructor
			return false;
		}
	}
	return true;
}

void IO::WriteSolutions(FILE *f)
{
	///////////////////////////Writing of the nodes////////////////////////////////
	fprintf(f, "\nSolutionSteps\t%d\n", db.number_solutions);
	for (int i = 0; i < db.number_solutions; i++)
		db.solution[i]->Write(f);			//Data writing 
}

void IO::WriteNodes(FILE *f)
{
	///////////////////////////Writing of the nodes////////////////////////////////
	fprintf(f, "\nNodes\t%d\n", db.number_nodes);
	for (int i = 0; i < db.number_nodes; i++)
		db.nodes[i]->Write(f);			//Data writing 
}

void IO::WritePoints(FILE *f)
{
	///////////////////////////Writing of the Points//////////////////////////////
	fprintf(f, "\nPoints\t%d\n", db.number_points);
	for (int i = 0; i < db.number_points; i++)
		db.points[i]->Write(f);			//Data writing 
}

void IO::WriteArcs(FILE *f)
{
	///////////////////////////Writing of the Arcs//////////////////////////////
	fprintf(f, "\nArcs\t%d\n", db.number_arcs);
	for (int i = 0; i < db.number_arcs; i++)
		db.arcs[i]->Write(f);			//Data writing 
}

void IO::WriteElements(FILE *f)
{
	///////////////////////////Writing of the elements///////////////////////////
	fprintf(f, "\nElements\t%d\n", db.number_elements);
	for (int i = 0; i < db.number_elements; i++)
		db.elements[i]->Write(f);			//Data writing
}
void IO::WriteParticles(FILE *f)
{
	///////////////////////////Writing of the particles///////////////////////////
	fprintf(f, "\nParticles\t%d\n", db.number_particles);
	for (int i = 0; i < db.number_particles; i++)
		db.particles[i]->Write(f);			//Data writing
}
void IO::WriteInitialConditions(FILE *f)
{
	///////////////////////////Writing of the initial conditions///////////////////
	fprintf(f, "\nInitialConditions\t%d\n", db.number_IC);
	for (int i = 0; i < db.number_IC; i++)
		db.IC[i]->Write(f);			//Data writing
}
void IO::WriteMaterials(FILE *f)
{
	///////////////////////////Writing of the materials///////////////////////////
	fprintf(f, "\nMaterials\t%d\n", db.number_materials);
	for (int i = 0; i < db.number_materials; i++)
		db.materials[i]->Write(f);		//Data writing 
}
void IO::WriteSections(FILE *f)
{
	///////////////////////////Writing of the sections//////////////////////////////
	fprintf(f, "\nSections\t%d\n", db.number_sections);
	for (int i = 0; i < db.number_sections; i++)
		db.sections[i]->Write(f);			//Data writing
}
void IO::WritePipeSections(FILE *f)
{
	///////////////////////////Writing of the sections of tubos/////////////////////
	fprintf(f, "\nPipeSections\t%d\n", db.number_pipe_sections);
	for (int i = 0; i < db.number_pipe_sections; i++)
		db.pipe_sections[i]->Write(f);	//Data writing
}
void IO::WriteShellSections(FILE *f)
{
	///////////////////////////Writing of the sections of tubos/////////////////////
	fprintf(f, "\nShellSections\t%d\n", db.number_shell_sections);
	for (int i = 0; i < db.number_shell_sections; i++)
		db.shell_sections[i]->Write(f);	//Data writing
}
void IO::WriteCoordinateSystems(FILE *f)
{
	///////////////////////////Writing of the CS//////////////////////////////////
	fprintf(f, "\nCoordinateSystems\t%d\n", db.number_CS);
	for (int i = 0; i < db.number_CS; i++)
		db.CS[i]->Write(f);			//Data writing
}

void IO::WriteRigidBodyData(FILE *f)
{
	///////////////////////////Writing of the CS//////////////////////////////////
	fprintf(f, "\nRigidBodyData\t%d\n", db.number_RB_data);
	for (int i = 0; i < db.number_RB_data; i++)
		db.RB_data[i]->Write(f);			//Data writing
}
void IO::WriteEnvironment(FILE *f)
{
	///////////////////////////Writing Environment/////////////////////////////
	fprintf(f, "\nEnvironment\n");
	db.environment->Write(f);
}

void IO::WriteMonitor(FILE *f)
{
	///////////////////////////Writing ////////////////////////////////
	fprintf(f, "\n");
	db.monitor->Write(f);			//Data writing 
}

void IO::WriteSolverOptions(FILE *f)
{
	///////////////////////////Writing SolverOptions////////////////////////
	fprintf(f, "\n");
	db.solver_options->Write(f);
}

void IO::WriteAnalyticalSurfaces(FILE *f)
{
	///////////////////////////Writing ////////////////////////////////
	fprintf(f, "\nAnalyticalSurfaces\t%d\n", db.number_analytical_surfaces);
	for (int i = 0; i < db.number_analytical_surfaces; i++)
		db.analytical_surfaces[i]->Write(f);			//Data writing 
}

void IO::WriteSurfaces(FILE *f)
{
	///////////////////////////Writing ////////////////////////////////
	fprintf(f, "\nSurfaces\t%d\n", db.number_surfaces);
	for (int i = 0; i < db.number_surfaces; i++)
		db.surfaces[i]->Write(f);			//Data writing 
}

void IO::WriteSplines(FILE *f)
{
	///////////////////////////Writing Splines////////////////////////
	fprintf(f, "\nSplines\t%d\n", db.number_splines);
	for (int i = 0; i < db.number_splines; i++)
		db.splines[i]->Write(f);
}

void IO::WriteLineRegions(FILE *f)
{
	///////////////////////////Writing ////////////////////////////////
	fprintf(f, "\nLineRegions\t%d\n", db.number_line_regions);
	for (int i = 0; i < db.number_line_regions; i++)
		db.line_regions[i]->Write(f);			//Data writing 
}

void IO::WriteSurfaceRegions(FILE *f)
{
	///////////////////////////Writing ////////////////////////////////
	fprintf(f, "\nSurfaceRegions\t%d\n", db.number_surface_regions);
	for (int i = 0; i < db.number_surface_regions; i++)
		db.surface_regions[i]->Write(f);			//Data writing 
}

void IO::WriteContacts(FILE *f)
{
	///////////////////////////Writing ////////////////////////////////
	fprintf(f, "\nContacts\t%d\n", db.number_contacts);
	for (int i = 0; i < db.number_contacts; i++)
		db.contacts[i]->Write(f);			//Data writing 
}

void IO::WriteNodeSets(FILE *f)
{
	///////////////////////////Writing NodeSets////////////////////////
	fprintf(f, "\nNodeSets\t%d\n", db.number_node_sets);
	for (int i = 0; i < db.number_node_sets; i++)
		db.node_sets[i]->Write(f);			//Data writing 
}

void IO::WriteSuperNodeSets(FILE *f)
{
	///////////////////////////Writing SuperNodeSets////////////////////////
	fprintf(f, "\nSuperNodeSets\t%d\n", db.number_super_node_sets);
	for (int i = 0; i < db.number_super_node_sets; i++)
		db.super_node_sets[i]->Write(f);			//Data writing 
}

void IO::WriteSurfaceSets(FILE *f)
{
	///////////////////////////Writing SurfaceSets////////////////////////
	fprintf(f, "\nSurfaceSets\t%d\n", db.number_surface_sets);
	for (int i = 0; i < db.number_surface_sets; i++)
		db.surface_sets[i]->Write(f);			//Data writing 
}

void IO::WriteElementSets(FILE *f)
{
	///////////////////////////Writing SurfaceSets////////////////////////
	fprintf(f, "\nElementSets\t%d\n", db.number_element_sets);
	for (int i = 0; i < db.number_element_sets; i++)
		db.element_sets[i]->Write(f);			//Data writing 
}

void IO::WriteLoads(FILE *f)
{
	///////////////////////////Writing of the Loads///////////////////////////
	fprintf(f, "\nLoads\t%d\n", db.number_loads);
	for (int i = 0; i < db.number_loads; i++)
		db.loads[i]->Write(f);		//Data writing 
}

void IO::WriteDisplacements(FILE *f)
{
	///////////////////////////Writing of the Displacements///////////////////////////
	fprintf(f, "\nDisplacements\t%d\n", db.number_displacements);
	for (int i = 0; i < db.number_displacements; i++)
		db.displacements[i]->Write(f);		//Data writing 
}

void IO::WriteConstraints(FILE *f)
{
	///////////////////////////Writing Constraints///////////////////////////
	fprintf(f, "\nConstraints\t%d\n", db.number_constraints);
	for (int i = 0; i < db.number_constraints; i++)
		db.constraints[i]->Write(f);	//Data writing
}

void IO::WriteSpecialConstraints(FILE *f)
{
	///////////////////////////Writing Special Constraints///////////////////////////
	fprintf(f, "\nSpecialConstraints\t%d\n", db.number_special_constraints);
	for (int i = 0; i < db.number_special_constraints; i++)
		db.special_constraints[i]->Write(f);	//Data writing
}

void IO::WriteSectionDetails(FILE *f)
{
	///////////////////////////Writing of the sections of tubos/////////////////////
	fprintf(f, "\nSectionDetails\t%d\n", db.number_section_details);
	for (int i = 0; i < db.number_section_details; i++)
		db.section_details[i]->Write(f);	//Data writing
}

void IO::WriteAerodynamicData(FILE *f)
{
	///////////////////////////Writing of the aerodynamicdata////////////////////
	fprintf(f, "\nAerodynamicData\t%d\n", db.number_aerodynamicdata);
	for (int i = 0; i < db.number_aerodynamicdata; i++)
		db.aerodynamic_data[i]->Write(f);			//Data writing 
}

void IO::WriteCADData(FILE *f)
{
	///////////////////////////Writing of the CADs///////////////////////////
	fprintf(f, "\nCADData\t%d\n", db.number_cad_data);
	for (int i = 0; i < db.number_cad_data; i++)
		db.cad_data[i]->Write(f);			//Data writing
}

void IO::WriteContactInterfaces(FILE *f)
{
	///////////////////////////Writing of the ContactInterfaces///////////////////////////
	fprintf(f, "\nContactInterfaces\t%d\n", db.number_contactinterfaces);
	for (int i = 0; i < db.number_contactinterfaces; i++)
		db.contactinterfaces[i]->Write(f);			//Data writing
}

void IO::WriteGeometries(FILE *f)
{
	///////////////////////////Writing of the Geometries///////////////////////////
	fprintf(f, "\nGeometries\t%d\n", db.number_geometries);
	for (int i = 0; i < db.number_geometries; i++)
		db.geometries[i]->Write(f);				//Data writing
}

void IO::WriteBodyGeometries(FILE *f)
{
	///////////////////////////Writing of the BodyGeometries///////////////////////////
	fprintf(f, "\nBodyGeometries\t%d\n", db.number_body_geometries);
	for (int i = 0; i < db.number_body_geometries; i++)
		db.body_geometries[i]->Write(f);		//Data writing
}

void IO::WriteBoundaries(FILE *f)
{
	///////////////////////////Writing of the Boundaries///////////////////////////
	fprintf(f, "\nBoundaries\t%d\n", db.number_boundaries);
	for (int i = 0; i < db.number_boundaries; i++)
		db.boundaries[i]->Write(f);			//Data writing
}

void IO::WriteBEMData(FILE *f)
{
	//////////////////////////Writing BEM////////////////////////////
	db.bem->Write(f);
}

void IO::WriteGeneralContactSearch(FILE *f)
{
	//////////////////////////Writing GCS////////////////////////////
	db.gcs->Write(f);
}

void IO::WriteConfigurationSave(FILE *f)
{
	//////////////////////////Writing ParticlePack////////////////////////////
	db.config_save->Write(f);
}

void IO::WriteConcomitantSolution(FILE *f)
{
	//////////////////////////Writing ConcomitantSolution///////////////////
	db.concomitant_solution->Write(f);
}

void IO::WriteConvergenceCriteria(FILE *f)
{
	//////////////////////////Writing Conv Criteria////////////////////////////
	fprintf(f, "\nConvergenceCriteria\n");
	db.conv_criteria->Write(f);
}

void IO::WritePostFiles(FILE *f)
{
	///////////////////////////Writing PostFiles////////////////////////
	fprintf(f, "\n");
	db.post_files->Write(f);
}


void IO::WritePSYCoupling(FILE *f)
{
	///////////////////////////Writing PSYCoupling////////////////////////
	fprintf(f, "\n");
	db.psy_coupling->Write(f);
}

void IO::WriteExecutionData(FILE* f)
{
	///////////////////////////Writing ExecutionData////////////////////////
	fprintf(f, "\n");
	db.execution_data->Write(f);
}

void IO::WriteSuperNodes(FILE *f)
{
	///////////////////////////Writing of the super nodes////////////////////////////////
	fprintf(f, "\nSuperNodes\t%d\n", db.number_super_nodes);
	for (int i = 0; i < db.number_super_nodes; i++)
		db.super_nodes[i]->Write(f);			//Data writing 
}
