#include "Load.h"

#include "Table.h"
#include "MathCode.h"

Load::Load()
{
}


Load::~Load()
{
}

double Load::GetValueAt(double t, int position)
{
	if (table != NULL)
		return table->GetValueAt(t, position);
	if (mcode != NULL)
		return mcode->GetValueAt(t, position);
	return 0.0;
}