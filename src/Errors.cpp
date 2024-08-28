#include "Errors.h"

#include "ConcomitantSolution.h"
#include "Solution.h"
#include "Node.h"
#include "Point.h"
#include "ArcCirc.h"
#include "Element.h"
#include "Particle.h"
#include "Boundary.h"
#include "BodyGeometry.h"
#include "Geometry.h"
#include "InitialCondition.h"
#include "Material.h"
#include "Section.h"
#include "PipeSection.h"
#include "ShellSection.h"
#include "CoordinateSystem.h"
#include "RigidBodyData.h"
#include "Monitor.h"
#include "Surface.h"
#include "AnalyticalSurface.h"
#include "LineRegion.h"
#include "SurfaceRegion.h"
#include "Contact.h"
#include "NodeSet.h"
#include "SuperNodeSet.h"
#include "SurfaceSet.h"
#include "ElementSet.h"
#include "Load.h"
#include "Displacement.h"
#include "Constraint.h"
#include "SpecialConstraint.h"
#include "AerodynamicData.h"
#include "CADData.h"
#include "ContactInterface.h"
#include "BEM.h"
#include "GeneralContactSearch.h"
#include "SuperNode.h"
#include "Spline.h"
#include "SectionDetails.h"
#include "Database.h"
//Variáveis globais
extern
Database db;

Errors::Errors()
{
}

Errors::~Errors()
{
}

bool Errors::CheckErrors()
{
	if (db.solution_exist == true)
		if (!CheckSolutions()) return false;
	if (db.nodes_exist == true)
		if (!CheckNodes()) return false;
	if (db.points_exist == true)
		if (!CheckPoints()) return false;
	if (db.arcs_exist == true)
		if (!CheckArcs()) return false;
	if (db.elements_exist == true)
		if (!CheckElements()) return false;
	if (db.particles_exist)
		if (!CheckParticles()) return false;
	if (db.boundaries_exist)
		if (!CheckBoundaries()) return false;
	if (db.IC_exist == true)
		if (!CheckInitialConditions()) return false;
	if (db.materials_exist == true)
		if (!CheckMaterials()) return false;
	if (db.sections_exist == true)
		if (!CheckSections()) return false;
	if (db.pipe_sections_exist == true)
		if (!CheckPipeSections()) return false;
	if (db.shell_sections_exist == true)
		if (!CheckShellSections()) return false;
	if (db.CS_exist == true)
		if (!CheckCoordinateSystems()) return false;
	if (db.RB_data_exist == true)
		if (!CheckRigidBodyData()) return false;
	if (db.environment_exist == true)
		if (!CheckEnvironment()) return false;
	if (db.monitor_exist == true)
		if (!CheckMonitor()) return false;
	if (db.solver_options_exist == true)
		if (!CheckSolverOptions()) return false;
	if (db.analytical_surfaces_exist == true)
		if (!CheckAnalyticalSurfaces()) return false;
	if (db.surfaces_exist == true)
		if (!CheckSurfaces()) return false;
	if (db.line_regions_exist == true)
		if (!CheckLineRegions()) return false;
	if (db.surface_regions_exist == true)
		if (!CheckSurfaceRegions()) return false;
	if (db.contacts_exist == true)
		if (!CheckContacts()) return false;
	if (db.node_sets_exist == true)
		if (!CheckNodeSets()) return false;
	if (db.super_node_sets_exist == true)
		if (!CheckSuperNodeSets()) return false;
	if (db.surface_sets_exist == true)
		if (!CheckSurfaceSets()) return false;
	if (db.element_sets_exist == true)
		if (!CheckElementSets()) return false;
	if (db.loads_exist == true)
		if (!CheckLoads()) return false;
	if (db.displacements_exist == true)
		if (!CheckDisplacements()) return false;
	if (db.constraints_exist == true)
		if (!CheckConstraints()) return false;
	if (db.special_constraints_exist == true)
		if (!CheckSpecialConstraints()) return false;
	if (db.section_details_exist == true)
		if (!CheckSectionDetails()) return false;
	if (db.aerodynamic_data_exist == true)
		if (!CheckAerodynamicData()) return false;
	if (db.cad_data_exist == true)
		if (!CheckCADData()) return false;
	if (db.contactinterfaces_exist == true)
		if (!CheckContacInterfaces()) return false;
	if (db.bem_exist == true)
		if (!CheckBEMData()) return false;
	if (db.concomitant_solution_exist == true)
		if (!CheckConcomitantSolution()) return false;
	if (db.gcs_exist == true)
		if (!CheckGeneralContactSearch()) return false;
	if (!CheckConvergenceCriteria()) return false;
	if (!CheckPostFiles()) return false;
	if (db.super_nodes_exist == true)
		if (!CheckSuperNodes()) return false;
	if (db.splines_exist == true)
		if (!CheckSplines()) return false;
	if (db.body_geometries_exist == true)
		if (!CheckBodyGeometries()) return false;
	if (db.geometries_exist == true)
		if (!CheckGometries()) return false;
	if (!CheckExecutionData()) return false;
	
	return true;
}

bool Errors::CheckSolutions()
{
	//Numbering
	for (int i = 0; i < db.number_solutions; i++)
	{
		if (db.solution[i]->solution_number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Solution data");
			return false;
		}
	}

	//Verifying if time sequence is plausible
	for (int i = 0; i < db.number_solutions-1; i++)
	{
		if (db.solution[i + 1]->end_time < db.solution[i]->end_time)
		{
			db.myprintf("Inconsistency in end time of Solution %d. Please, check the input file.\n", i+2);
			return false;
		}	
	} 
	return true;
}

bool Errors::CheckNodes()
{
	for (int i = 0; i < db.number_nodes; i++)
	{
		if (db.nodes[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Nodes data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckPoints()
{
	for (int i = 0; i < db.number_points; i++)
	{
		if (db.points[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Points data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckArcs()
{
	for (int i = 0; i < db.number_arcs; i++)
	{
		if (db.arcs[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Arcs data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckElements()
{
	for (int i = 0; i < db.number_elements; i++)
	{
		if (db.elements[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Elements data");
			return false;
		}
		if (!db.elements[i]->Check())
		{
			db.myprintf("Inconsistency in Element %d. Please, check the input file.\n", i+1);
			return false;
		}
			
	}
	return true;
}

bool Errors::CheckParticles()
{
	for (int i = 0; i < db.number_particles; i++)
	{
		if (db.particles[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Particles data");
			return false;
		}
		if (!db.particles[i]->Check())
		{
			db.myprintf("Inconsistency in Particle %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}


bool Errors::CheckBoundaries()
{
	for (int i = 0; i < db.number_boundaries; i++)
	{
		if (db.boundaries[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Boundaries data");
			return false;
		}
		if (!db.boundaries[i]->Check())
		{
			db.myprintf("Inconsistency in Boundary %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckBodyGeometries()
{
	for (int i = 0; i < db.number_body_geometries; i++)
	{
		if (db.body_geometries[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Body geometries data");
			return false;
		}
		if (!db.body_geometries[i]->Check())
		{
			db.myprintf("Inconsistency in Body Geometry %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckGometries()
{
	for (int i = 0; i < db.number_geometries; i++)
	{
		if (db.geometries[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Geometries data");
			return false;
		}
		if (!db.geometries[i]->Check())
		{
			db.myprintf("Inconsistency in Geometry %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckInitialConditions()
{
	for (int i = 0; i < db.number_IC; i++)
	{
		if (db.IC[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Initial Conditions data");
			return false;
		}
		if (!db.IC[i]->Check())
		{
			db.myprintf("Inconsistency in Initial Condition %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckMaterials()
{
	for (int i = 0; i < db.number_materials; i++)
	{
		if (db.materials[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Materials data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckSections()
{
	for (int i = 0; i < db.number_sections; i++)
	{
		if (db.sections[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Sections data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckPipeSections()
{
	for (int i = 0; i < db.number_pipe_sections; i++)
	{
		if (db.pipe_sections[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Pipe Sections data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckShellSections()
{
	for (int i = 0; i < db.number_shell_sections; i++)
	{
		if (db.shell_sections[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Shell Sections data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckCoordinateSystems()
{
	for (int i = 0; i < db.number_CS; i++)
	{
		if (db.CS[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Coordinate Systems data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckRigidBodyData()
{
	for (int i = 0; i < db.number_RB_data; i++)
	{
		if (db.RB_data[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Rigid Body Data data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckEnvironment()
{
	return true;
}

bool Errors::CheckMonitor()
{
	for (int i = 0; i < db.monitor->nodes.size(); i++)
	{
		if (db.monitor->nodes[i] > db.number_nodes)
		{
			db.myprintf("Error in %s. Please, check the input file.\n", "Monitor, in Nodes requested for monitoring");
			return false;
		}
	}
	for (int i = 0; i < db.monitor->elements.size(); i++)
	{
		if (db.monitor->elements[i] > db.number_elements)
		{
			db.myprintf("Error in %s. Please, check the input file.\n", "Monitor, in Elements requested for monitoring");
			return false;
		}
	}
	for (int i = 0; i < db.monitor->contacts.size(); i++)
	{
		if (db.monitor->contacts[i] > db.number_contacts)
		{
			db.myprintf("Error in %s. Please, check the input file.\n", "Monitor, in Contacts requested for monitoring");
			return false;
		}
	}
	for (int i = 0; i < db.monitor->node_sets.size(); i++)
	{
		if (db.monitor->node_sets[i] > db.number_node_sets)
		{
			db.myprintf("Error in %s. Please, check the input file.\n", "Monitor, in Node Sets requested for monitoring");
			return false;
		}
	}
	return true;
}

bool Errors::CheckSolverOptions()
{
	return true;
}

bool Errors::CheckAnalyticalSurfaces()
{
	for (int i = 0; i < db.number_analytical_surfaces; i++)
	{
		if (db.analytical_surfaces[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Analytical Surfaces data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckSurfaces()
{
	for (int i = 0; i < db.number_surfaces; i++)
	{
		if (db.surfaces[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Surfaces data");
			return false;
		}
		if (!db.surfaces[i]->Check())
		{
			db.myprintf("Inconsistency in Surface %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckLineRegions()
{
	for (int i = 0; i < db.number_line_regions; i++)
	{
		if (db.line_regions[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Line Regions data");
			return false;
		}
		for (int e = 0; e < db.line_regions[i]->n_elements; e++)
		{
			if (db.line_regions[i]->elements[e] > db.number_elements)
			{
				db.myprintf("Error in %s. Please, check the input file.\n", "Line Regions data");
				return false;
			}
		}
	}
	return true;
}

bool Errors::CheckSurfaceRegions()
{
	for (int i = 0; i < db.number_surface_regions; i++)
	{
		if (db.surface_regions[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Surface Regions data");
			return false;
		}
		for (int e = 0; e < db.surface_regions[i]->n_elements; e++)
		{
			if (db.surface_regions[i]->elements[e] > db.number_elements)
			{
				db.myprintf("Error in %s. Please, check the input file.\n", "Surface Regions data");
				return false;
			}
		}
	}
	return true;
}

bool Errors::CheckContacts()
{
	for (int i = 0; i < db.number_contacts; i++)
	{
		if (db.contacts[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Contacts data");
			return false;
		}
		if (!db.contacts[i]->Check())
		{
			db.myprintf("Inconsistency in Contact %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckNodeSets()
{
	for (int i = 0; i < db.number_node_sets; i++)
	{
		if (db.node_sets[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Node Sets data");
			return false;
		}
		for (int e = 0; e < db.node_sets[i]->n_nodes; e++)
		{
			if (db.node_sets[i]->node_list[e] > db.number_nodes)
			{
				db.myprintf("Error in %s. Please, check the input file.\n", "Node Sets data");
				return false;
			}
		}
	}
	return true;
}

bool Errors::CheckSuperNodeSets()
{
	for (int i = 0; i < db.number_super_node_sets; i++)
	{
		if (db.super_node_sets[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Super Node Sets data");
			return false;
		}
		for (int e = 0; e < db.super_node_sets[i]->n_super_nodes; e++)
		{
			if (db.super_node_sets[i]->super_node_list[e] > db.number_super_nodes)
			{
				db.myprintf("Error in %s. Please, check the input file.\n", "Super Node Sets data");
				return false;
			}
		}
	}
	return true;
}

bool Errors::CheckSurfaceSets()
{
	for (int i = 0; i < db.number_surface_sets; i++)
	{
		if (db.surface_sets[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Surface Sets data");
			return false;
		}
		for (int e = 0; e < db.surface_sets[i]->n_surf; e++)
		{
			if (db.surface_sets[i]->surf_list[e] > db.number_surfaces)
			{
				db.myprintf("Error in %s. Please, check the input file.\n", "Surface Sets data");
				return false;
			}
		}
	}
	return true;
}

bool Errors::CheckElementSets()
{
	for (int i = 0; i < db.number_element_sets; i++)
	{
		if (db.element_sets[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Element Sets data");
			return false;
		}
		for (int e = 0; e < db.element_sets[i]->n_el; e++)
		{
			if (db.element_sets[i]->el_list[e] > db.number_elements)
			{
				db.myprintf("Error in %s. Please, check the input file.\n", "Element Sets data");
				return false;
			}
		}
	}
	return true;
}

bool Errors::CheckLoads()
{
	for (int i = 0; i < db.number_loads; i++)
	{
		if (db.loads[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Loads data");
			return false;
		}
		if (!db.loads[i]->Check())
		{
			db.myprintf("Inconsistency in Load %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckDisplacements()
{
	for (int i = 0; i < db.number_displacements; i++)
	{
		if (db.displacements[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Displacements data");
			return false;
		}
		if (!db.displacements[i]->Check())
		{
			db.myprintf("Inconsistency in Displacement %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckConstraints()
{
	for (int i = 0; i < db.number_constraints; i++)
	{
		if (db.constraints[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Constraints data");
			return false;
		}
		if (!db.constraints[i]->Check())
		{
			db.myprintf("Inconsistency in Constraint %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckSpecialConstraints()
{
	for (int i = 0; i < db.number_special_constraints; i++)
	{
		if (db.special_constraints[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Special Constraints data");
			return false;
		}
		if (!db.special_constraints[i]->Check())
		{
			db.myprintf("Inconsistency in Special Constraint %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckSectionDetails()
{
	for (int i = 0; i < db.number_section_details; i++)
	{
		if (db.section_details[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Section Details data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckAerodynamicData()
{
	for (int i = 0; i < db.number_aerodynamicdata; i++)
	{
		if (db.aerodynamic_data[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Aerodynamic data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckCADData()
{
	for (int i = 0; i < db.number_cad_data; i++)
	{
		if (db.cad_data[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "CAD data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckContacInterfaces()
{
	for (int i = 0; i < db.number_contactinterfaces; i++)
	{
		if (db.contactinterfaces[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "ContacInterfaces data");
			return false;
		}
		if (!db.contactinterfaces[i]->Check())
		{
			db.myprintf("Inconsistency in ContacInterface %d. Please, check the input file.\n", i + 1);
			return false;
		}
	}
	return true;
}

bool Errors::CheckBEMData()
{
	if (!db.bem->Check())
	{
		db.myprintf("Inconsistency in BEM data.\n");
		return false;
	}
	return true;
}

bool Errors::CheckConcomitantSolution()
{
	if (!db.concomitant_solution->Check())
	{
		db.myprintf("Inconsistency in ConcomitantSolution data.\n");
		return false;
	}
	return true;
}


bool Errors::CheckConvergenceCriteria()
{
	return true;
}

bool Errors::CheckGeneralContactSearch()
{
	if (!db.gcs->Check())
	{
		db.myprintf("Missing data for GeneralContactSearch evaluation. Please, check contact interfaces.\n");
		return false;
	}
	return true;
	
}

bool Errors::CheckPostFiles()
{
	return true;
}

bool Errors::CheckSuperNodes()
{
	for (int i = 0; i < db.number_super_nodes; i++)
	{
		if (db.super_nodes[i]->ID - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Super Nodes data");
			return false;
		}
	}
	return true;
}

bool Errors::CheckExecutionData()
{
	return true;
}

bool Errors::CheckSplines()
{
	for (int i = 0; i < db.number_splines; i++)
	{
		if (db.splines[i]->number - 1 != i)
		{
			db.myprintf("Inconsistency in %s numbering. Please, check the input file.\n", "Splines data");
			return false;
		}
		if (db.splines[i]->nodeset > db.number_node_sets)
		{
			db.myprintf("Inconsistency in %s definition. Please, check the input file.\n", "Spline nodeset");
			return false;
		}
		if ((db.node_sets[db.splines[i]->nodeset - 1]->n_nodes - 1) % 2 != 0) 
		{
			db.myprintf("Inconsistency in %s definition for spline. Please, check the input file.\n", "number of nodes in Nodeset");
			return false;
		}
	}
	return true;
}
