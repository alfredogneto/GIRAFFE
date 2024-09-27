#include "Tetrahedron.h"


Tetrahedron::Tetrahedron()
{
	CAD_ID = 0;
	verticesIDs[0] = 0;
	verticesIDs[1] = 0;
	verticesIDs[2] = 0;
	verticesIDs[3] = 0;
	edgesIDs[0] = 0;
	edgesIDs[1] = 0;
	edgesIDs[2] = 0;
	edgesIDs[3] = 0;
	edgesIDs[4] = 0;
	edgesIDs[5] = 0;
	ID = 0;												
}

Tetrahedron::~Tetrahedron()
{
}
void Tetrahedron::Print(FILE *f)
{
	fprintf(f, "Tetrahedron %d\t%d\t%d\t%d\t%d\n", ID - 1, verticesIDs[0] - 1, verticesIDs[1] - 1, verticesIDs[2] - 1, verticesIDs[3] - 1);
}