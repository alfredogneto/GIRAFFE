#include "SecUserDefined.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

SecUserDefined::SecUserDefined()
{
	//Variables
	GA = 0.0;
	EA = 0.0;
	ES1 = 0.0;
	ES2 = 0.0;
	EI11 = 0.0;
	EI22 = 0.0;
	EI12 = 0.0;
	GS1 = 0.0;
	GS2 = 0.0;
	GS1S = 0.0;
	GS2S = 0.0;
	GJT = 0.0;
	J11 = 0.0;
	J22 = 0.0;
	J12 = 0.0;
	A = 0.0;
	SC = Matrix(2);
	BC = Matrix(2);
	Rho = 0.0;
	sec_details_ID = 0;
	sec_details = NULL;

	aerodynamicdataID = 0;
	AC = Matrix(2);
	aero_length = 0;
}

SecUserDefined::~SecUserDefined()
{
	if (sec_details != NULL)
		delete sec_details;
}

bool SecUserDefined::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "GA"))
	{
		fscanf(f, "%s", s);
		GA = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "EA"))
	{
		fscanf(f, "%s", s);
		EA = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "ES1"))
	{
		fscanf(f, "%s", s);
		ES1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "ES2"))
	{
		fscanf(f, "%s", s);
		ES2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "EI11"))
	{
		fscanf(f, "%s", s);
		EI11 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "EI22"))
	{
		fscanf(f, "%s", s);
		EI22 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "EI12"))
	{
		fscanf(f, "%s", s);
		EI12 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "GS1"))
	{
		fscanf(f, "%s", s);
		GS1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "GS2"))
	{
		fscanf(f, "%s", s);
		GS2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "GS1S"))
	{
		fscanf(f, "%s", s);
		GS1S = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "GS2S"))
	{
		fscanf(f, "%s", s);
		GS2S = atof(s);
	}
	else
		return false;
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "GJT"))
	{
		fscanf(f, "%s", s);
		GJT = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "J11"))
	{
		fscanf(f, "%s", s);
		J11 = atof(s);
	}
	else
		return false;
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "J22"))
	{
		fscanf(f, "%s", s);
		J22 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "J12"))
	{
		fscanf(f, "%s", s);
		J12 = atof(s);
	}
	else
		return false;
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "A"))
	{
		fscanf(f, "%s", s);
		A = atof(s);
	}
	else
		return false;
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "SC"))
	{
		fscanf(f, "%s", s);
		SC(0, 0) = atof(s);
		fscanf(f, "%s", s);
		SC(1, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "BC"))
	{
		fscanf(f, "%s", s);
		BC(0, 0) = atof(s);
		fscanf(f, "%s", s);
		BC(1, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Rho"))
	{
		fscanf(f, "%s", s);
		Rho = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "SD"))
	{
		fscanf(f, "%s", s);
		sec_details_ID = atoi(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "AD"))
	{
		fscanf(f, "%s", s);
		aerodynamicdataID = atoi(s);

		fscanf(f, "%s", s);
		if (!strcmp(s, "AC"))
		{
			fscanf(f, "%s", s);
			AC(0, 0) = atof(s);
			fscanf(f, "%s", s);
			AC(1, 0) = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "AeroLength"))
		{
			fscanf(f, "%s", s);
			aero_length = atof(s);
		}
		else
			return false;
	}
	else
		fsetpos(f, &pos);	//volta à posição anterior

	return true;
}

void SecUserDefined::Write(FILE *f)
{
	/*
	GA = 0.0;
	EA = 0.0;
	ES1 = 0.0;
	ES2 = 0.0;
	EI11 = 0.0;
	EI22 = 0.0;
	EI12 = 0.0;
	GS1 = 0.0;
	GS2 = 0.0;
	GS1S = 0.0;
	GS2S = 0.0;
	GJT = 0.0;
	J11 = 0.0;
	J22 = 0.0;
	J12 = 0.0;
	A = 0.0;
	SC = Matrix(2);
	BC = Matrix(2);
	Rho = 0.0;
	SD = 1;
	*/

	fprintf(f, "UserDefined\t%d\nGA\t%.6e\nEA\t%.6e\nES1\t%.6e\nES2\t%.6e\nEI11\t%.6e\nEI22\t%.6e\nEI12\t%.6e\nGS1\t%.6e\nGS2\t%.6e\nGS1S\t%.6e\nGS2S\t%.6e\nGJT\t%.6e\nJ11\t%.6e\nJ22\t%.6e\nJ12\t%.6e\nA\t%.6e\nSC\t%.6e\t%.6e\nBC\t%.6e\t%.6e\nRho\t%.6e\nSD\t%d\nAD\t%dAC\t%.6e\t%.6e\nAeroLength\t%.6e\n",
		number, GA, EA, ES1, ES2, EI11, EI22, EI12, GS1, GS2, GS1S, GS2S, GJT, J11, J22, J12, A, SC(0,0), SC(1,0), BC(0,0), BC(1,0), Rho,sec_details_ID,aerodynamicdataID, AC(0,0), AC(1,0), aero_length);
}
void SecUserDefined::PreCalc()
{
	//Atribuição dos section details à seção transversal
	if (sec_details_ID <= db.number_section_details)
	{
		//Verifica o tipo de SD atribuído e aloca de acordo
		if (typeid(*db.section_details[sec_details_ID - 1]) == typeid(SolidSection))
		{
			SolidSection* ptr_sd = static_cast<SolidSection*>(db.section_details[sec_details_ID - 1]);
			sec_details = new SolidSection(*ptr_sd);
		}
		if (typeid(*db.section_details[sec_details_ID - 1]) == typeid(MultiCellSection))
		{
			MultiCellSection* ptr_sd = static_cast<MultiCellSection*>(db.section_details[sec_details_ID - 1]);
			sec_details = new MultiCellSection(*ptr_sd);
		}
		
	}
	else
	{
		db.myprintf("Error. SD number associated to UserDefined cross section number %d does not exist. It was ignored.\n", number);
		sec_details = new SolidSection();
	}	
}