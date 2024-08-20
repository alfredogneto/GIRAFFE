#pragma once

class Errors
{
public:
	Errors();
	~Errors();
	bool CheckErrors();		//Checa erros de inconsistência do arquivo de entrada e alerta o usuário caso encontre

	//Funções individuais para checar erros em cada um dos ítens
	bool CheckSolutions();
	bool CheckNodes();
	bool CheckPoints();
	bool CheckArcs();
	bool CheckElements();
	bool CheckParticles();
	bool CheckBoundaries();
	bool CheckInitialConditions();
	bool CheckMaterials();
	bool CheckSections();
	bool CheckPipeSections();
	bool CheckShellSections();
	bool CheckCoordinateSystems();
	bool CheckRigidBodyData();
	bool CheckEnvironment();
	bool CheckMonitor();
	bool CheckSolverOptions();
	bool CheckAnalyticalSurfaces();
	bool CheckSurfaces();
	bool CheckLineRegions();
	bool CheckSurfaceRegions();
	bool CheckContacts();
	bool CheckNodeSets();
	bool CheckSuperNodeSets();
	bool CheckSurfaceSets();
	bool CheckElementSets();
	bool CheckLoads();
	bool CheckDisplacements();
	bool CheckConstraints();
	bool CheckSpecialConstraints();
	bool CheckSectionDetails();
	bool CheckAerodynamicData();
	bool CheckCADData();
	bool CheckContacInterfaces();
	bool CheckBEMData();
	bool CheckConcomitantSolution();
	bool CheckConvergenceCriteria();
	bool CheckGeneralContactSearch();
	bool CheckPostFiles();
	bool CheckSuperNodes();
	bool CheckExecutionData();
	bool CheckSplines();
	bool CheckBodyGeometries();
	bool CheckGometries();
};