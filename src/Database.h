#pragma once
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

//Elementos
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

//Materiais
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

//Contatos
#include "Contact.h"
#include "LRLR.h"
#include "GeneralPLR.h"
#include "NSSS.h"
#include "SSSS.h"
#include "SPSP.h"

//GeneralContactSearch
#include "GeneralContactSearch.h"

//Seções transversais
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

//Partículas
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

//Deslocamentos prescritos
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

//OPENMP
#include <omp.h>
#include <process.h>

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