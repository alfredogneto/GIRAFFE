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

//Elementos
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

//Materiais
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

//Contatos
class Contact;
class LRLR;
class GeneralPLR;
class NSSS;
class SSSS;
class SPSP;

//GeneralContactSearch
class GeneralContactSearch;

//Seções transversais
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

//Partículas
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

//Deslocamentos prescritos
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
	void PrintPtr(double** ptr, int lines, int columns);				//Imprime na tela o conteúdo do ponteiro ptr (matriz de duas dimensões)
	void PrintPtr(double* ptr, int lines);								//Imprime na tela o conteúdo do ponteiro ptr (matriz de uma dimensão)
	double EvaluateBoundingBoxDiag();	//Calcula uma dimensão característica da geometria do modelo - de todos os nós existentes - com base em suas posições atuais
	void PreCalc();						//Realiza pré-cálculo
	void EvaluateStartEndTimes();		//Raliza pré-cálculo de start times com base na sequencia de solutions

	//Variables and objects
	char file_name[1000];			//Nome do job
	char folder_name[1000];			//Pasta em que estão sendo salvos arquivos da simulação

	int number_solutions;			//Número de soluções
	Solution** solution;			//Vetor de soluções

	int number_nodes;				//Número de nós
	Node** nodes;					//Vetor de nós

	int number_points;				//Número de pontos
	Point** points;					//Vetor de pontos

	int number_arcs;				//Número de arcos
	ArcCirc** arcs;					//Vetor de arcos

	int number_elements;			//Número de elementos
	Element** elements;				//Vetor de elementos

	int number_particles;			//Número de partículas
	Particle** particles;			//Vetor de partículas

	int number_IC;					//Número de condições iniciais
	InitialCondition** IC;			//Vetor de condições iniciais

	int number_materials;			
	Material** materials;			//Vetor de materiais

	int number_sections;
	Section** sections;				//Vetor de sections

	int number_pipe_sections;
	PipeSection** pipe_sections;	//Vetor de pipe_sections

	int number_shell_sections;
	ShellSection** shell_sections;	//Vetor de pipe_sections

	int number_CS;
	CoordinateSystem** CS;			//Vetor de coordinate systems

	int number_RB_data;
	RigidBodyData** RB_data;		//Vetor de Rigid body data

	Environment* environment;				//Environment
	
	Monitor* monitor;						//Monitor de resultados
	
	SolverOptions* solver_options;			//SolverOptions

	int number_analytical_surfaces;
	AnalyticalSurface** analytical_surfaces;		//Vetor de superfícies analíticas

	int number_surfaces;
	Surface** surfaces;								//Vetor de superfícies

	int number_splines;								//Número de splines
	Spline** splines;								//Vetor de splines

	int number_line_regions;
	LineRegion** line_regions;						//Vetor de line regions

	int number_surface_regions;
	SurfaceRegion** surface_regions;				//Vetor de surface regions

	int number_contacts;
	Contact** contacts;								//Vetor de contatos

	int number_node_sets;
	NodeSet** node_sets;							//Vetor de node sets

	int number_super_node_sets;
	SuperNodeSet** super_node_sets;					//Vetor de super node sets

	int number_surface_sets;
	SurfaceSet** surface_sets;						//Vetor de surface sets

	int number_element_sets;
	ElementSet** element_sets;						//Vetor de element sets

	int number_loads;
	Load** loads;									//Vetor de loads

	int number_displacements;
	Displacement** displacements;					//Vetor de displacements

	int number_constraints;
	Constraint** constraints;						//Vetor de constraints

	int number_special_constraints;
	SpecialConstraint** special_constraints;		//Vetor de special constraints

	int number_section_details;
	SectionDetails** section_details;				//Vetor de section_details

	int number_aerodynamicdata;
	AerodynamicData** aerodynamic_data;				//Vetor de dados aerodinâmicos

	int number_cad_data;							//Número de CADs
	CADData** cad_data;								//Vetor de CADs

	int number_contactinterfaces;					//Número de interfaces de contato
	ContactInterface** contactinterfaces;			//Vetor de interfaces de contato

	int number_boundaries;
	Boundary** boundaries;							//Vetor de boundaries

	int number_super_nodes;							//Número de super nós
	SuperNode** super_nodes;						//Vetor de super nós

	int number_body_geometries;						//Número de body contact boundaries
	BodyGeometry** body_geometries;					//Vetor de body contact boundaries

	int number_geometries;							//Número de geometries
	Geometry** geometries;							//Vetor de geometries

	BEM* bem;										//Blade element momentum method

	GeneralContactSearch* gcs;						//General Contact Search

	ConfigurationSave* config_save;					//Particle pack routines

	ConcomitantSolution* concomitant_solution;		//Concomitant solution

	ConvergenceCriteria* conv_criteria;				//Convergence criteria
	PostFiles* post_files;							//Pós processamento (arquivos)
	ExecutionData* execution_data;					//Dados de execução do GIRAFFE

	//Variáveis booleanas de controle
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

	bool flag_nGL_changed;		//Flag - indica que houve mudança nos graus de liberdade
	int n_GL_free;				//GLs livres
	int n_GL_fixed;				//GLs prescritos
	int number_GLs_node;		//Número de GLs por nó

	//Variáveis para execução do programa
	double last_converged_time;				//Último instante convergido
	double current_time_step;				//Incremento de tempo atual
	int current_solution_number;			//Atual solution
	int current_iteration_number;			//Iteration
	bool plot_times;
	
	int size_AA;							//Estimativa de coef. não nulos da matriz AA
	int size_BB;							//Estimativa de coef. não nulos da matriz BB
	int size_AB;							//Estimativa de coef. não nulos da matriz AB

	SparseMatrix global_stiffness_AA;		//Matriz de rigidez global
	SparseMatrix global_stiffness_BB;		//Matriz de rigidez global
	SparseMatrix global_stiffness_AB;		//Matriz de rigidez global
	SparseMatrix global_stiffness_BA;		//Matriz de rigidez global

	SparseMatrix global_mass_AA;			//Matriz de massa global - utilizada somente quando há análise modal		
	SparseMatrix global_damping_AA;			//Matriz de amrotecimento global - utilizada somente quando há análise modal

	Matrix global_P_A;				//Esforço desbalanceado global (após a resolução do sistema linear, também armazena os deslocamentos generalizados)
	Matrix global_I_A;				//Esforços internos globais - gl livres
	Matrix global_P_B;				//Esforços internos globais - gl fixos
	Matrix global_X_B;				//Vetor de deslocamentos impostos nos GL do sistema

	//Matrix global_ABS_P_A;			//Vetor para avaliação do valor absoluto do "fluxo" prescrito em cada DOF
	//Matrix global_COUNT_ABS_P_A;	//Vetor que conta quantas contribuições foram colocadas em cada DOF
	//Matrix global_ABS_P_B;			//Vetor para avaliação do valor absoluto do "fluxo" prescrito em cada DOF
	//Matrix global_COUNT_ABS_P_B;	//Vetor que conta quantas contribuições foram colocadas em cada DOF

	//Info para saída de resultados de elementos - número de resultados por elemento
	int n_element_results;

	FILE* console_output;	//Stream para salvar o console_output para um arquivo de texto
	char version[10];		//Versão do GIRAFFE
	char name[1000];		//Nome do arquivo do output da solução

	//PSY
	PSYCoupling* psy_coupling;
	bool psy_coupling_exist;
};