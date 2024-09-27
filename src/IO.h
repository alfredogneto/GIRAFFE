#pragma once
#include "Database.h"

//Variaveis globais
extern 
Database db;
class IO
{
public:
	IO(void);
	~IO(void);

	bool ReadFile(int argc, char* argv[]);
	void WriteFile();

	bool ReadSolutions(FILE *f);
	bool ReadNodes(FILE *f);
	bool ReadPoints(FILE *f);
	bool ReadArcs(FILE *f);
	bool ReadElements(FILE *f);
	bool ReadParticles(FILE *f);
	bool ReadInitialConditions(FILE *f);
	bool ReadMaterials(FILE *f);
	bool ReadSections(FILE *f);
	bool ReadPipeSections(FILE *f);
	bool ReadShellSections(FILE *f);
	bool ReadCoordinateSystems(FILE *f);
	bool ReadRigidBodyData(FILE *f);
	bool ReadEnvironment(FILE *f);
	bool ReadMonitor(FILE *f);
	bool ReadSolverOptions(FILE *f);
	bool ReadAnalyticalSurfaces(FILE *f);
	bool ReadSurfaces(FILE *f);
	bool ReadSplines(FILE *f);
	bool ReadLineRegions(FILE *f);
	bool ReadSurfaceRegions(FILE *f);
	bool ReadContacts(FILE *f);
	bool ReadNodeSets(FILE *f);
	bool ReadSuperNodeSets(FILE *f);
	bool ReadSurfaceSets(FILE *f);
	bool ReadElementSets(FILE *f);
	bool ReadLoads(FILE *f);
	bool ReadDisplacements(FILE *f);
	bool ReadConstraints(FILE *f);
	bool ReadSpecialConstraints(FILE *f);
	bool ReadSectionDetails(FILE *f);
	bool ReadAerodynamicData(FILE *f);
	bool ReadCADData(FILE *f);
	bool ReadContactInterfaces(FILE *f);
	bool ReadBoundaries(FILE *f);
	bool ReadBEMData(FILE *f);
	bool ReadGeneralContactSearch(FILE *f);
	bool ReadConfigurationSave(FILE *f);
	bool ReadConcomitantSolution(FILE *f);
	bool ReadConvergenceCriteria(FILE *f);
	bool ReadPostFiles(FILE *f);
	bool ReadPSYCoupling(FILE *f);
	bool ReadExecutionData(FILE* f);
	bool ReadSuperNodes(FILE *f);
	bool ReadGeometries(FILE *f);
	bool ReadBodyGeometries(FILE *f);

	void WriteSolutions(FILE *f);
	void WriteNodes(FILE *f);
	void WritePoints(FILE *f);
	void WriteArcs(FILE *f);
	void WriteElements(FILE *f);
	void WriteParticles(FILE *f);
	void WriteInitialConditions(FILE *f);
	void WriteMaterials(FILE *f);
	void WriteSections(FILE *f);
	void WritePipeSections(FILE *f);
	void WriteShellSections(FILE *f);
	void WriteCoordinateSystems(FILE *f);
	void WriteRigidBodyData(FILE *f);
	void WriteEnvironment(FILE *f);
	void WriteMonitor(FILE *f);
	void WriteSolverOptions(FILE *f);
	void WriteAnalyticalSurfaces(FILE *f);
	void WriteSurfaces(FILE *f);
	void WriteSplines(FILE* f);
	void WriteLineRegions(FILE *f);
	void WriteSurfaceRegions(FILE *f);
	void WriteContacts(FILE *f);
	void WriteNodeSets(FILE *f);
	void WriteSuperNodeSets(FILE *f);
	void WriteSurfaceSets(FILE *f);
	void WriteElementSets(FILE *f);
	void WriteLoads(FILE *f);
	void WriteDisplacements(FILE *f);
	void WriteConstraints(FILE *f);
	void WriteSpecialConstraints(FILE *f);
	void WriteSectionDetails(FILE *f);
	void WriteAerodynamicData(FILE *f);
	void WriteCADData(FILE *f);
	void WriteContactInterfaces(FILE *f);
	void WriteBoundaries(FILE *f);
	void WriteBEMData(FILE *f);
	void WriteGeneralContactSearch(FILE *f);
	void WriteConfigurationSave(FILE *f);
	void WriteConcomitantSolution(FILE *f);
	void WriteConvergenceCriteria(FILE *f);
	void WritePostFiles(FILE *f);
	void WritePSYCoupling(FILE* f);
	void WriteExecutionData(FILE* f);
	void WriteSuperNodes(FILE *f);
	void WriteGeometries(FILE *f);
	void WriteBodyGeometries(FILE *f);
};

bool ReadComment(FILE *f, char* s, int dim_char);	//Lê comentarios - retorna o stream no ponto após leitura de comentario
void TryComment(FILE *f);							//Tenta ler comentarios. Retorna o stream no ponto em que a próxima leitura não e um comentario
void WriteDOFTable(FILE *f);						//Escreve tabela com DOFs