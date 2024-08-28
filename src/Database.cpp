#include "Database.h"
#include <stdarg.h>

#include "ConvergenceCriteria.h"
#include "PostFiles.h"
#include "SolverOptions.h"
#include "ExecutionData.h"
#include "Node.h"
#include "CADData.h"
#include "BodyGeometry.h"
#include "Load.h"
#include "Constraint.h"
#include "Displacement.h"
#include "SpecialConstraint.h"
#include "ArcCirc.h"
#include "Section.h"
#include "Element.h"
#include "AnalyticalSurface.h"
#include "Surface.h"
#include "Spline.h"
#include "Contact.h"
#include "Particle.h"
#include "Boundary.h"
#include "PipeSection.h"
#include "ContactInterface.h"
#include "PSYCoupling.h"
#include "GeneralContactSearch.h"
#include "Monitor.h"
#include "Modal.h"
#include "SplineElement.h"

Database::Database()
{
	//Maximum number of DOFs/node
	number_GLs_node = 6;
	//Versão do GIRAFFE
	sprintf(version,"2.0.135");
	
	number_solutions = 0;
	number_nodes = 0;
	number_super_nodes = 0;
	number_points = 0;
	number_arcs = 0;
	number_elements = 0;
	number_particles = 0;
	number_IC = 0;
	number_materials = 0;
	number_sections = 0;
	number_pipe_sections = 0;
	number_shell_sections = 0;
	number_CS = 0;
	number_RB_data = 0;
	number_analytical_surfaces = 0;
	number_surfaces = 0;
	number_splines = 0;
	number_line_regions = 0;
	number_surface_regions = 0;
	number_contacts = 0;
	number_node_sets = 0;
	number_super_node_sets = 0;
	number_surface_sets = 0;
	number_element_sets = 0;
	number_loads = 0;
	number_displacements = 0;
	number_constraints = 0;
	number_special_constraints = 0;
	number_section_details = 0;
	number_aerodynamicdata = 0;
	number_cad_data = 0;
	number_contactinterfaces = 0;
	number_boundaries = 0;
	number_super_nodes = 0;
	number_body_geometries = 0;
	number_geometries = 0;

	solution = NULL;
	nodes = NULL;
	super_nodes = NULL;
	points = NULL;
	arcs = NULL;
	elements = NULL;
	particles = NULL;
	IC = NULL;
	materials = NULL;
	sections = NULL;
	pipe_sections = NULL;
	shell_sections = NULL;
	CS = NULL;
	RB_data = NULL;
	environment = NULL;
	monitor = NULL;
	analytical_surfaces = NULL;
	surfaces = NULL;
	splines = NULL;
	line_regions = NULL;
	surface_regions = NULL;
	contacts = NULL;
	node_sets = NULL;
	super_node_sets = NULL;
	surface_sets = NULL;
	element_sets = NULL;
	loads = NULL;
	displacements = NULL;
	constraints = NULL;
	special_constraints = NULL;
	section_details = NULL;
	aerodynamic_data = NULL;
	cad_data = NULL;
	contactinterfaces = NULL;
	boundaries = NULL;
	super_nodes = NULL;
	body_geometries = NULL;
	geometries = NULL;

	bem = NULL;
	gcs = NULL;
	config_save = NULL;
	concomitant_solution = NULL;
	psy_coupling = NULL;

	//Objetos que sempre existem no Database
	conv_criteria = new ConvergenceCriteria();
	post_files = new PostFiles();
	solver_options = new SolverOptions();
	execution_data = new ExecutionData();

	//Variaveis booleanas de controle
	solution_exist = false;
	nodes_exist = false;
	super_nodes_exist = false;
	points_exist = false;
	arcs_exist = false;
	elements_exist = false;
	particles_exist = false;
	IC_exist = false;
	materials_exist = false;
	sections_exist = false;
	pipe_sections_exist = false;
	shell_sections_exist = false;
	CS_exist = false;
	RB_data_exist = false;
	environment_exist = false;
	monitor_exist = false;
	solver_options_exist = true;
	analytical_surfaces_exist = false;
	surfaces_exist = false;
	splines = false;
	line_regions_exist = false;
	surface_regions_exist = false;
	contacts_exist = false;
	node_sets_exist = false;
	super_node_sets_exist = false;
	surface_sets_exist = false;
	element_sets_exist = false;
	loads_exist = false;
	displacements_exist = false;
	constraints_exist = false;
	special_constraints_exist = false;
	section_details_exist = false;
	aerodynamic_data_exist = false;
	cad_data_exist = false;
	contactinterfaces_exist = false;
	boundaries_exist = false;
	super_nodes_exist = false;
	body_geometries_exist = false;
	geometries_exist = false;

	bem_exist = false;
	gcs_exist = false;
	config_save_exist = false;
	concomitant_solution_exist = false;
	psy_coupling_exist = false;

	flag_nGL_changed = true;
	n_GL_free = 0;
	n_GL_fixed = 0;
	
	last_converged_time = 0.0;				//ultimo instante convergido
	current_time_step = 0.0;				//Incremento de tempo atual
	current_solution_number = 0;			//Atual solution
	current_iteration_number = 0;
	plot_times = false; 
	
	size_AA = 0;
	size_BB = 0;
	size_AB = 0;

	n_element_results = 15;
	console_output = NULL;
}

Database::~Database()
{
	if (solution != NULL)
	{
		for (int i = 0; i < this->number_solutions; i++)
			delete solution[i];
		delete[] solution;
	}

	if (nodes != NULL)
	{
		for (int i = 0; i < this->number_nodes; i++)
			delete nodes[i];
		delete[] nodes;
	}

	if (points != NULL)
	{
		for (int i = 0; i < this->number_points; i++)
			delete points[i];
		delete[] points;
	}

	if (arcs != NULL)
	{
		for (int i = 0; i < this->number_arcs; i++)
			delete arcs[i];
		delete[] arcs;
	}
	
	if (elements != NULL)
	{
		for (int i = 0; i < this->number_elements; i++)
			delete elements[i];
		delete[] elements;
	}

	if (particles != NULL)
	{
		for (int i = 0; i < this->number_particles; i++)
			delete particles[i];
		delete[] particles;
	}

	if (IC != NULL)
	{
		for (int i = 0; i < this->number_IC; i++)
			delete IC[i];
		delete[] IC;
	}

	if (materials != NULL)
	{
		for (int i = 0; i < this->number_materials; i++)
			delete materials[i];
		delete[] materials;
	}
	
	if (sections != NULL)
	{
		for (int i = 0; i < this->number_sections; i++)
			delete sections[i];
		delete[] sections;
	}
	
	if (pipe_sections != NULL)
	{
		for (int i = 0; i < this->number_pipe_sections; i++)
			delete pipe_sections[i];
		delete[] pipe_sections;
	}

	if (shell_sections != NULL)
	{
		for (int i = 0; i < this->number_shell_sections; i++)
			delete shell_sections[i];
		delete[] shell_sections;
	}
	
	if (CS != NULL)
	{
		for (int i = 0; i < this->number_CS; i++)
			delete CS[i];
		delete[] CS;
	}
	
	if (RB_data != NULL)
	{
		for (int i = 0; i < this->number_RB_data; i++)
			delete RB_data[i];
		delete[] RB_data;
	}

	if (environment != NULL)
		delete environment;

	if (monitor != NULL)
		delete monitor;

	if (solver_options != NULL)
		delete solver_options;

	if (analytical_surfaces != NULL)
	{
		for (int i = 0; i < this->number_analytical_surfaces; i++)
			delete analytical_surfaces[i];
		delete[] analytical_surfaces;
	}

	if (surfaces != NULL)
	{
		for (int i = 0; i < this->number_surfaces; i++)
			delete surfaces[i];
		delete[] surfaces;
	}

	if (splines != NULL)
	{
		for (int i = 0; i < this->number_splines; i++)
			delete splines[i];
		delete[] splines;
	}

	if (line_regions != NULL)
	{
		for (int i = 0; i < this->number_line_regions; i++)
			delete line_regions[i];
		delete[] line_regions;
	}

	if (surface_regions != NULL)
	{
		for (int i = 0; i < this->number_surface_regions; i++)
			delete surface_regions[i];
		delete[] surface_regions;
	}

	if (contacts != NULL)
	{
		for (int i = 0; i < this->number_contacts; i++)
			delete contacts[i];
		delete[] contacts;
	}

	if (node_sets != NULL)
	{
		for (int i = 0; i < this->number_node_sets; i++)
			delete node_sets[i];
		delete[] node_sets;
	}

	if (super_node_sets != NULL)
	{
		for (int i = 0; i < this->number_super_node_sets; i++)
			delete super_node_sets[i];
		delete[] super_node_sets;
	}

	if (surface_sets != NULL)
	{
		for (int i = 0; i < this->number_surface_sets; i++)
			delete surface_sets[i];
		delete[] surface_sets;
	}

	if (element_sets != NULL)
	{
		for (int i = 0; i < this->number_element_sets; i++)
			delete element_sets[i];
		delete[] element_sets;
	}

	if (loads != NULL)
	{
		for (int i = 0; i < this->number_loads; i++)
			delete loads[i];
		delete[] loads;
	}

	if (displacements != NULL)
	{
		for (int i = 0; i < this->number_displacements; i++)
			delete displacements[i];
		delete[] displacements;
	}

	if (constraints != NULL)
	{
		for (int i = 0; i < this->number_constraints; i++)
			delete constraints[i];
		delete[] constraints;
	}

	if (special_constraints != NULL)
	{
		for (int i = 0; i < this->number_special_constraints; i++)
			delete special_constraints[i];
		delete[] special_constraints;
	}

	if (section_details != NULL)
	{
		for (int i = 0; i < this->number_section_details; i++)
			delete section_details[i];
		delete[] section_details;
	}

	if (aerodynamic_data != NULL)
	{
		for (int i = 0; i < this->number_aerodynamicdata; i++)
			delete aerodynamic_data[i];
		delete[] aerodynamic_data;
	}

	if (cad_data != NULL)
	{
		for (int i = 0; i < this->number_cad_data; i++)
			delete cad_data[i];
		delete[] cad_data;
	}

	if (contactinterfaces != NULL)
	{
		for (int i = 0; i < this->number_contactinterfaces; i++)
			delete contactinterfaces[i];
		delete[] contactinterfaces;
	}

	if (boundaries != NULL)
	{
		for (int i = 0; i < this->number_boundaries; i++)
			delete boundaries[i];
		delete[] boundaries;
	}

	if (super_nodes != NULL)
	{
		for (int i = 0; i < this->number_super_nodes; i++)
			delete super_nodes[i];
		delete[] super_nodes;
	}

	if (body_geometries != NULL)
	{
		for (int i = 0; i < this->number_body_geometries; i++)
			delete body_geometries[i];
		delete[] body_geometries;
	}

	if (geometries != NULL)
	{
		for (int i = 0; i < this->number_geometries; i++)
			delete geometries[i];
		delete[] geometries;
	}

	if (bem != NULL)
		delete bem;

	if (gcs != NULL)
		delete gcs;

	if (config_save != NULL)
		delete config_save;

	if (concomitant_solution != NULL)
		delete concomitant_solution;

	if (psy_coupling != NULL)
		delete psy_coupling;

	delete conv_criteria;
	delete post_files;
	delete execution_data;
}

int Database::myprintf(const char* format, ...)
{
	int ret_status = 0;
	va_list args;

	if (execution_data->print_console)
	{
		va_start(args, format);
		ret_status = vprintf(format, args);		//impressão na tela
		va_end(args);
	}
	
	if (console_output != NULL && execution_data->print_file)
	{
		va_start(args, format);
		ret_status = vfprintf(console_output, format, args);		//impressão no arquivo de texto
		va_end(args);
	}
	return ret_status;
}

//Imprime na tela o conteudo do ponteiro ptr (matriz de duas dimensões)
void Database::PrintPtr(double* ptr, int lines)
{
	myprintf("\n");
	for (int i = 0; i < lines; i++)
	{
		if (abs(ptr[i]) > 1e-8)
			myprintf("%.6e ", ptr[i]);
		else
			myprintf("%.2f ", ptr[i]);
		myprintf("\n");
	}
	myprintf("\n");
}
//Imprime na tela o conteudo do ponteiro ptr (matriz de duas dimensões)
void Database::PrintPtr(double** ptr, int lines, int columns)
{
	myprintf("\n");
	for (int i = 0; i < lines; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			if (abs(ptr[i][j]) > 1e-8)
				myprintf("%.6e ", ptr[i][j]);
			else
				myprintf("%.2f ", ptr[i][j]);
		}
		myprintf("\n");
	}
	myprintf("\n");
}

//Calcula uma dimensão caracteristica da geometria do modelo - de todos os nós existentes - com base em suas posições atuais
double Database::EvaluateBoundingBoxDiag()
{
	//Verify bounding box size
	double max_x = 0;
	double min_x = 0;
	double max_y = 0;
	double min_y = 0;
	double max_z = 0;
	double min_z = 0;
	//Determinação do bounding box
	for (int i = 0; i < number_nodes; i++)
	{
		if (nodes[i]->copy_coordinates[0]>max_x)
			max_x = nodes[i]->copy_coordinates[0];
		if (nodes[i]->copy_coordinates[0]<min_x)
			min_x = nodes[i]->copy_coordinates[0];
		if (nodes[i]->copy_coordinates[1]>max_y)
			max_y = nodes[i]->copy_coordinates[1];
		if (nodes[i]->copy_coordinates[1]<min_y)
			min_y = nodes[i]->copy_coordinates[1];
		if (nodes[i]->copy_coordinates[2]>max_z)
			max_z = nodes[i]->copy_coordinates[2];
		if (nodes[i]->copy_coordinates[2]<min_z)
			min_z = nodes[i]->copy_coordinates[2];
	}
	//calcula a diagonal do bounding box e retorna
	return sqrt((max_x - min_x)*(max_x - min_x) + (max_y - min_y)*(max_y - min_y) + (max_z - min_z)*(max_z - min_z));
}
//Realiza pre-calculo em todos os elementos
void Database::PreCalc()
{
#pragma omp for
	for (int i = 0; i < number_cad_data; i++)
		cad_data[i]->PreCalc();
	for (int i = 0; i < number_arcs; i++)
		arcs[i]->PreCalc();
	for (int i = 0; i < number_sections; i++)
		sections[i]->PreCalc();
	for (int i = 0; i < number_elements; i++)
		elements[i]->PreCalc();
	for (int i = 0; i < number_analytical_surfaces; i++)
		analytical_surfaces[i]->PreCalc();
	for (int i = 0; i < number_surfaces; i++)
		surfaces[i]->PreCalc();
	for (int i = 0; i < number_splines; i++)
		splines[i]->PreCalc();
	for (int i = 0; i < number_splines; i++)
		for (int j = 0; j < splines[i]->size_sp_elements; j++)
			splines[i]->sp_element[j]->PreCalc();
	for (int i = 0; i < number_contacts; i++)
		contacts[i]->PreCalc();
	for (int i = 0; i < number_particles; i++)
		particles[i]->PreCalc();
	for (int i = 0; i < number_boundaries; i++)
		boundaries[i]->PreCalc();
	//for (int i = 0; i < number_geometries; i++)		//Called inside body_geometries[i]->PreCalc();
	//	geometries[i]->PreCalc();
	for (int i = 0; i < number_body_geometries; i++)
		body_geometries[i]->PreCalc();
	for (int i = 0; i < number_loads; i++)
		loads[i]->PreCalc();
	for (int i = 0; i < number_constraints; i++)
		constraints[i]->PreCalc();
	for (int i = 0; i < number_displacements; i++)
		displacements[i]->PreCalc();
	for (int i = 0; i < number_special_constraints; i++)
		special_constraints[i]->PreCalc();
	for (int i = 0; i < number_pipe_sections; i++)
		pipe_sections[i]->PreCalc();
	for (int i = 0; i < number_contactinterfaces; i++)
		contactinterfaces[i]->PreCalc();

	solver_options->PreCalc();

	//PSY
	if (psy_coupling_exist)
		psy_coupling->PreCalc();
	if (gcs_exist)
		gcs->PreCalc();

	//Aloca soluções do PostFiles
	post_files->AllocFiles(number_solutions + 1);	//the last is only the pvd file pointing to all other solutions
	if (monitor_exist == true)
		monitor->AllocFiles();
}
void Database::EvaluateStartEndTimes()
{
	double current_start = 0.0;
	for (int i = 0; i < number_solutions; i++)
	{
		solution[i]->start_time = current_start;
		if (typeid(*solution[i]) == typeid(Modal))
			solution[i]->end_time = current_start;
		current_start = solution[i]->end_time;
	}
}
