#pragma once
#include "Matrix.h"
#include "SparseMatrix.h"

class Node;
class SuperNode;
class PipeSection;
class RigidBodyData;
class CoordinateSystem;
class Environment;

class Monitor;
class PostFiles;
class SolverOptions;

class NodeSet;
class SurfaceSet;
class ElementSet;
class SuperNodeSet;

class InitialCondition;
class ConvergenceCriteria;

class AnalyticalSurface;
class Plane;
class LineRegion;
class SurfaceRegion;

//Elements
class Element;
class Beam_1;
class Pipe_1;
class Shell_1;
class Solid_1;
class SpringDashpot_1;
class Mass_1;
class RigidBody_1;
class Truss_1;
class TwoNodeConnector_1;

//Materials
class Material;
class Hooke;
class ElasticPlasticIsoHardening;
class Orthotropic;

//Point & Curves
class Point;
class ArcCirc;

//Surfaces
class Surface;
class RigidTriangularSurface_1;
class RigidOscillatorySurface_1;
class FlexibleTriangularSurface_2;
class FlexibleSECylinder_1;
class FlexibleArcExtrusion_1;
class RigidArcRevolution_1;
class RigidNURBS_1;

//Geometries
class Geometry;
class SECylinder;
class ArcExtrusion;
class ArcRevolution;

//BodyGeometry
class BodyGeometry;

//Splines
class Spline;
class SplineElement;

//Contacts
class Contact;
class LRLR;
class GeneralPLR;
class NSSS;
class SSSS;
class SPSP;

//GeneralContactSearch
class GeneralContactSearch;

//Sections transversais
class Section;
class SecGeneral;
class SecRectangle;
class SecSuperEllipse;
class SecTube;
class SecUserDefined;
class SecHelicalFiber;

//Shell sections
class ShellSection;
class ShellSectionHomogeneous;
class ShellSectionComposite;

//Special constraints
class SpecialConstraint;
class SameDisplacement;
class Hinge;
class UniversalJoint;
class SameRotation;
class RigidNodeSet;
class TranslationalJoint;
class NodalConstraintDOF;

//Particles
class Particle;
class Sphere;
class Polyhedron;
class NURBSParticle;
class VEMPolyhedron;

//Boundaries
class Boundary;
class STLBoundary;

//Surface pairs (for general contact search)
class SurfacePairGeneralContact;
class RigidTriangularFace_RigidTriangularFace;
class FlexibleTriangularFace_FlexibleTriangularFace;
class FlexibleTriangularFace_RigidTriangularFace;

//Section details
class SectionDetails;
class SolidSection;
class MultiCellSection;

//AerodynamicData
class AerodynamicData;
class BEM;

//Solutions
class Solution;
class Static;
class Dynamic;
class Modal;
class ConcomitantSolution;
class ExplicitDynamic;

//Loads
class Load;
class NodalLoad;
class NodalFollowerLoad;
class ShellLoad;
class PipeLoad;
class SuperNodalLoad;

//Displacements prescritos
class Displacement;
class NodalDisplacement;
class DisplacementField;

//Constraints
class Constraint;
class NodalConstraint;
class SuperNodalConstraint;

//PSY
class PSYCoupling;

//CADData
class CADData;
class STLSurface;
class NURBSSurface;

//ContactInterfaces
class ContactInterface;
class Interface_0;
class Interface_1;

//BoundingVolumes
class BoundingVolume;
class BoundingSphere;
class BoundingCylinder;
class BoundingTriangularBox;
class BoundingBoxAxesAligned;
class BoundingBoxAxesOriented;

class Encoding;
class ExecutionData;
class ConfigurationSave;



//OPENMP
//#include <omp.h>
//#include <process.h>

class Database
{
public:
	Database();
	~Database();

	//Methods
	int myprintf(const char* format, ...);
	void PrintPtr(double** ptr, int lines, int columns);				//Prints in the tela the content of the pointer ptr (matrix of two dimensions)
	void PrintPtr(double* ptr, int lines);								//Prints in the tela the content of the pointer ptr (matrix of a dimension)
	double EvaluateBoundingBoxDiag();	//Calculates a dimension caracteristica of the geometry of the modelo - of all the nodes existentes - with base in its positions current
	void PreCalc();						//Performs pre-calculation
	void EvaluateStartEndTimes();		//Raliza pre-calculation of start times with base in the sequencia of solutions

	//Variables and objects
	char file_name[1000];			//Name of the job
	char folder_name[1000];			//Folder in that estão being saved files of the simulation

	int number_solutions;			//Number of solutions
	Solution** solution;			//Vector of solutions

	int number_nodes;				//Number of nodes
	Node** nodes;					//Vector of nodes

	int number_points;				//Number of points
	Point** points;					//Vector of points

	int number_arcs;				//Number of arcos
	ArcCirc** arcs;					//Vector of arcos

	int number_elements;			//Number of elements
	Element** elements;				//Vector of elements

	int number_particles;			//Number of particles
	Particle** particles;			//Vector of particles

	int number_IC;					//Number of initial conditions
	InitialCondition** IC;			//Vector of initial conditions

	int number_materials;			
	Material** materials;			//Vector of materials

	int number_sections;
	Section** sections;				//Vector of sections

	int number_pipe_sections;
	PipeSection** pipe_sections;	//Vector of pipe_sections

	int number_shell_sections;
	ShellSection** shell_sections;	//Vector of pipe_sections

	int number_CS;
	CoordinateSystem** CS;			//Vector of coordinate systems

	int number_RB_data;
	RigidBodyData** RB_data;		//Vector of Rigid body data

	Environment* environment;				//Environment
	
	Monitor* monitor;						//Monitor of results
	
	SolverOptions* solver_options;			//SolverOptions

	int number_analytical_surfaces;
	AnalyticalSurface** analytical_surfaces;		//Vector of surfaces analiticas

	int number_surfaces;
	Surface** surfaces;								//Vector of surfaces

	int number_splines;								//Number of splines
	Spline** splines;								//Vector of splines

	int number_line_regions;
	LineRegion** line_regions;						//Vector of line regions

	int number_surface_regions;
	SurfaceRegion** surface_regions;				//Vector of surface regions

	int number_contacts;
	Contact** contacts;								//Vector of contacts

	int number_node_sets;
	NodeSet** node_sets;							//Vector of node sets

	int number_super_node_sets;
	SuperNodeSet** super_node_sets;					//Vector of super node sets

	int number_surface_sets;
	SurfaceSet** surface_sets;						//Vector of surface sets

	int number_element_sets;
	ElementSet** element_sets;						//Vector of element sets

	int number_loads;
	Load** loads;									//Vector of loads

	int number_displacements;
	Displacement** displacements;					//Vector of displacements

	int number_constraints;
	Constraint** constraints;						//Vector of constraints

	int number_special_constraints;
	SpecialConstraint** special_constraints;		//Vector of special constraints

	int number_section_details;
	SectionDetails** section_details;				//Vector of section_details

	int number_aerodynamicdata;
	AerodynamicData** aerodynamic_data;				//Vector of data aerodinamicos

	int number_cad_data;							//Number of CADs
	CADData** cad_data;								//Vector of CADs

	int number_contactinterfaces;					//Number of interfaces of contact
	ContactInterface** contactinterfaces;			//Vector of interfaces of contact

	int number_boundaries;
	Boundary** boundaries;							//Vector of boundaries

	int number_super_nodes;							//Number of super nodes
	SuperNode** super_nodes;						//Vector of super nodes

	int number_body_geometries;						//Number of body contact boundaries
	BodyGeometry** body_geometries;					//Vector of body contact boundaries

	int number_geometries;							//Number of geometries
	Geometry** geometries;							//Vector of geometries

	BEM* bem;										//Blade element momentum method

	GeneralContactSearch* gcs;						//General Contact Search

	ConfigurationSave* config_save;					//Particle pack routines

	ConcomitantSolution* concomitant_solution;		//Concomitant solution

	ConvergenceCriteria* conv_criteria;				//Convergence criteria
	PostFiles* post_files;							//Pós processamento (files)
	ExecutionData* execution_data;					//Data of execution of the GIRAFFE

	//Variables booleanas of control
	bool solution_exist;
	bool nodes_exist;
	bool super_nodes_exist;
	bool points_exist;
	bool arcs_exist;
	bool elements_exist;
	bool particles_exist;
	bool IC_exist;
	bool materials_exist;
	bool sections_exist;
	bool pipe_sections_exist;
	bool shell_sections_exist;
	bool CS_exist;
	bool RB_data_exist;
	bool environment_exist;
	bool monitor_exist;
	bool solver_options_exist;
	bool analytical_surfaces_exist;
	bool surfaces_exist;
	bool splines_exist;
	bool line_regions_exist;
	bool surface_regions_exist;
	bool contacts_exist;
	bool node_sets_exist;
	bool super_node_sets_exist;
	bool surface_sets_exist;
	bool element_sets_exist;
	bool loads_exist;
	bool displacements_exist;
	bool constraints_exist;
	bool special_constraints_exist;
	bool section_details_exist;
	bool aerodynamic_data_exist;
	bool cad_data_exist;
	bool contactinterfaces_exist;
	bool boundaries_exist;
	bool body_geometries_exist;
	bool geometries_exist;

	bool bem_exist;
	bool gcs_exist;
	bool config_save_exist;
	bool concomitant_solution_exist;

	bool flag_nGL_changed;		//Flag - indicates that there was change in the degrees of freedom
	int n_GL_free;				//GLs free
	int n_GL_fixed;				//GLs prescritos
	int number_GLs_node;		//Number of GLs by node

	//Variables for execution of the programa
	double last_converged_time;				//last instant converged
	double current_time_step;				//Time increment current
	int current_solution_number;			//Current solution
	int current_iteration_number;			//Iteration
	bool plot_times;
	
	int size_AA;							//Estimativa of coef. not null of the matrix AA
	int size_BB;							//Estimativa of coef. not null of the matrix BB
	int size_AB;							//Estimativa of coef. not null of the matrix AB

	SparseMatrix global_stiffness_AA;		//Stiffness matrix global
	SparseMatrix global_stiffness_BB;		//Stiffness matrix global
	SparseMatrix global_stiffness_AB;		//Stiffness matrix global
	SparseMatrix global_stiffness_BA;		//Stiffness matrix global

	SparseMatrix global_mass_AA;			//Mass matrix global - used only when there is analysis modal		
	SparseMatrix global_damping_AA;			//Matrix of amrotecimento global - used only when there is analysis modal

	Matrix global_P_A;				//Internal force unbalanced global (after the resolução of the linear system, tambem armazena the displacements generalizados)
	Matrix global_I_A;				//Internal forces internal global - gl free
	Matrix global_P_B;				//Internal forces internal global - gl fixed
	Matrix global_X_B;				//Displacement vector impostos in the GL of the system

	//Matrix global_ABS_P_A;			//Vector for evaluation of the value absoluto of the "fluxo" prescrito in each DOF
	//Matrix global_COUNT_ABS_P_A;	//Vector that conta quantas contributions were colocadas in each DOF
	//Matrix global_ABS_P_B;			//Vector for evaluation of the value absoluto of the "fluxo" prescrito in each DOF
	//Matrix global_COUNT_ABS_P_B;	//Vector that conta quantas contributions were colocadas in each DOF

	//Info for output of results of elements - number of results by element
	int n_element_results;

	FILE* console_output;	//Stream for save the console_output for a file of text
	char version[10];		//Version of the GIRAFFE
	char name[1000];		//File name of the output of the solution

	//PSY
	PSYCoupling* psy_coupling;
	bool psy_coupling_exist;
};