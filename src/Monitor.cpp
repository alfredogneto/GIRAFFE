#include "Monitor.h"
#include <direct.h>

#include "Node.h"
#include "Element.h"
#include "Contact.h"
#include "NodeSet.h"
#include "ElementSet.h"
#include "Particle.h"
#include "SpecialConstraint.h"
#include "GeneralContactSearch.h"
#include "ConvergenceCriteria.h"
#include "Pipe_1.h"

#include"Database.h"
#include "IO.h"
//Variaveis globais
extern
Database db;


Monitor::Monitor()
{
	nodes.clear();
	elements.clear();
	contacts.clear();
	node_sets.clear();
	global.clear();
	particles.clear();
	special_constraints.clear();

	f_nodes = NULL;
	f_elements = NULL;
	f_contacts = NULL;
	f_node_sets = NULL;
	f_global = NULL;
	f_particles = NULL;
	f_special_constraints = NULL;
	
	first_record = true;
	sample = 1;

	alloced_files = false;

	contact_special_output = false;
	print_times = true;

	monitor_nodes_exists = false;
	monitor_elements_exists = false;
	monitor_contacts_exists = false;
	monitor_node_sets_exists = false;
	monitor_global_exists = false;
	monitor_particles_exists = false;
	monitor_special_constraints_exists = false;
	monitor_userdef_exists = false;
	monitor_elemuserdef_exists = false;

	userdef_entries.clear();
	userdef_elem_entries.clear();
}

Monitor::~Monitor()
{
	FlushFiles();

	nodes.clear();
	elements.clear();
	contacts.clear();
	node_sets.clear();
	global.clear();
	particles.clear();
	special_constraints.clear();
	userdef_entries.clear();
	userdef_elem_entries.clear();
}

void Monitor::AllocFiles()
{
	FlushFiles();
	f_nodes = new FILE*[(int)nodes.size()];
	f_elements = new FILE*[(int)elements.size()];
	f_contacts = new FILE*[(int)contacts.size()];
	f_node_sets = new FILE*[(int)node_sets.size()];
	f_particles = new FILE*[(int)particles.size()];
	f_global = new FILE*[1];
	f_special_constraints = new FILE*[(int)special_constraints.size()];
	alloced_files = true;
}

void Monitor::FlushFiles()
{
	if (alloced_files == true)
	{
		if (f_nodes != NULL)
			delete[]f_nodes;
		if (f_elements != NULL)
			delete[]f_elements;
		if (f_contacts != NULL)
			delete[]f_contacts;
		if (f_node_sets != NULL)
			delete[]f_node_sets;
		if (f_particles != NULL)
			delete[]f_particles;
		if (f_global != NULL)
			delete[]f_global;
		if (f_special_constraints != NULL)
			delete[]f_special_constraints;
	
		f_nodes = NULL;
		f_elements = NULL;
		f_contacts = NULL;
		f_node_sets = NULL;
		f_particles = NULL;
		f_global = NULL;
		f_special_constraints = NULL;
		
		alloced_files = false;
	}
}

bool Monitor::Read(FILE *f)
{
	char s[1000];
	//Verifica a palavra chave "Sample"
	fscanf(f, "%s", s);
	if (!strcmp(s, "Sample"))
	{
		fscanf(f, "%s", s);
		sample = atoi(s);
	}
	else
		return false;

	bool option_OK = false;
	bool flag_continue = true;
	fpos_t pos;
	while (flag_continue == true)
	{
		option_OK = false;
		TryComment(f);
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF)
			return true;
		if (!strcmp(s, "MonitorNodes"))
		{
			if (!ReadIntTable(f,&nodes))
				return false;
			option_OK = true;
			monitor_nodes_exists = true;
		}
		if (!strcmp(s, "MonitorElements"))
		{
			if (!ReadIntTable(f, &elements))
				return false;
			option_OK = true;
			monitor_elements_exists = true;
		}
		if (!strcmp(s, "MonitorContacts"))
		{
			if (!ReadIntTable(f, &contacts))
				return false;
			option_OK = true;
			monitor_contacts_exists = true;
		}
		if (!strcmp(s, "MonitorContactsSpecial"))
		{
			contact_special_output = true;
			if (!ReadIntTable(f, &contacts))
				return false;
			option_OK = true;
			monitor_contacts_exists = true;
		}
		if (!strcmp(s, "MonitorNodeSets"))
		{
			if (!ReadIntTable(f, &node_sets))
				return false;
			option_OK = true;
			monitor_node_sets_exists = true;
		}
		if (!strcmp(s, "MonitorParticles"))
		{
			if (!ReadIntTable(f, &particles))
				return false;
			option_OK = true;
			monitor_particles_exists = true;
		}
		if (!strcmp(s, "MonitorGlobal"))
		{
			option_OK = true;
			monitor_global_exists = true;
		}

		if (!strcmp(s, "MonitorSpecialConstraints"))
		{
			if (!ReadIntTable(f, &special_constraints))
				return false;
			option_OK = true;
			monitor_special_constraints_exists = true;
		}
		
		if (!strcmp(s, "MonitorNodesUserDefined"))
		{
			if (!ReadUserDefEntry(f))
				return false;
			option_OK = true;
			monitor_userdef_exists = true;
		}

		if (!strcmp(s, "MonitorElementsUserDefined"))
		{
			if (!ReadElemUserDefEntry(f))
				return false;
			option_OK = true;
			monitor_elemuserdef_exists = true;
		}

		if (option_OK == true)
			flag_continue = true;
		else
		{
			fsetpos(f, &pos);
			flag_continue = false;
		}
	}
	return true;
}

bool Monitor::ReadIntTable(FILE *f, vector<int> *data_table)
{
	char s[1000];
	fpos_t pos;
	fgetpos(f, &pos);
	bool flag_not_digit = false;
	fscanf(f, "%s", s);
	while (flag_not_digit == false)
	{
		data_table->push_back((atoi(s)));
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF)
			return true;
		if (!isdigit(s[0]))
			flag_not_digit = true;
	}
	fsetpos(f, &pos);
	return true;
}

bool Monitor::ReadUserDefEntry(FILE *f)
{
	char s[1000];
	UserDefMonitorEntry entry;

	// Ler "Nodes" ou "NodeSet"
	if (fscanf(f, "%s", s) == EOF) return false;

	if (!strcmp(s, "Nodes"))
	{
		entry.use_node_set = false;
	}
	else if (!strcmp(s, "NodeSet"))
	{
		entry.use_node_set = true;
	}
	else
	{
		return false;
	}

	// Ler IDs inteiros at� encontrar token n�o-d�gito
	fpos_t pos;
	fgetpos(f, &pos);
	if (fscanf(f, "%s", s) == EOF) return false;
	while (isdigit(s[0]))
	{
		entry.ids.push_back(atoi(s));
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF) return false;
	}

	// Token n�o-d�gito deve ser "Parameters"
	if (strcmp(s, "Parameters") != 0) return false;

	// Ler par�metros
	if (!ReadUserDefParams(f, entry.params)) return false;

	// Verificar se há TimeStart opcional
	fpos_t pos2;
	char s2[1000];
	fgetpos(f, &pos2);
	if (fscanf(f, "%s", s2) != EOF)
	{
		if (!strcmp(s2, "TimeStart"))
		{
			if (fscanf(f, "%s", s2) == EOF) return false;
			entry.time_start = atof(s2);
		}
		else
		{
			fsetpos(f, &pos2);
		}
	}

	if (entry.ids.empty()) return false;

	userdef_entries.push_back(entry);
	return true;
}

bool Monitor::ReadUserDefParams(FILE *f, UserDefMonitorParams& p)
{
	char s[1000];
	fpos_t pos;
	bool found = false;

	while (true)
	{
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF) break;

		bool is_param = true;
		if      (!strcmp(s, "X"))   p.X   = true;
		else if (!strcmp(s, "Y"))   p.Y   = true;
		else if (!strcmp(s, "Z"))   p.Z   = true;
		else if (!strcmp(s, "dX"))  p.dX  = true;
		else if (!strcmp(s, "dY"))  p.dY  = true;
		else if (!strcmp(s, "dZ"))  p.dZ  = true;
		else if (!strcmp(s, "ddX")) p.ddX = true;
		else if (!strcmp(s, "ddY")) p.ddY = true;
		else if (!strcmp(s, "ddZ")) p.ddZ = true;
		else if (!strcmp(s, "FX"))  p.FX  = true;
		else if (!strcmp(s, "FY"))  p.FY  = true;
		else if (!strcmp(s, "FZ"))  p.FZ  = true;
		else if (!strcmp(s, "MX"))  p.MX  = true;
		else if (!strcmp(s, "MY"))  p.MY  = true;
		else if (!strcmp(s, "MZ"))  p.MZ  = true;
		else is_param = false;

		if (!is_param)
		{
			fsetpos(f, &pos);
			break;
		}
		found = true;
	}
	return found;
}

bool Monitor::ReadElemUserDefEntry(FILE *f)
{
	char s[1000];
	UserDefElemMonitorEntry entry;

	// Ler "Elements" ou "ElementSet"
	if (fscanf(f, "%s", s) == EOF) return false;

	if (!strcmp(s, "Elements"))
	{
		entry.use_element_set = false;
	}
	else if (!strcmp(s, "ElementSet"))
	{
		entry.use_element_set = true;
	}
	else
	{
		return false;
	}

	// Ler IDs inteiros at� encontrar token n�o-d�gito
	fpos_t pos;
	fgetpos(f, &pos);
	if (fscanf(f, "%s", s) == EOF) return false;
	while (isdigit(s[0]))
	{
		entry.ids.push_back(atoi(s));
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF) return false;
	}

	// Token n�o-d�gito deve ser "Parameters"
	if (strcmp(s, "Parameters") != 0) return false;

	// Ler par�metros
	if (!ReadElemUserDefParams(f, entry.params)) return false;

	// Verificar se h� TimeStart opcional
	fpos_t pos2;
	fgetpos(f, &pos2);
	char s2[1000];
	if (fscanf(f, "%s", s2) != EOF)
	{
		if (!strcmp(s2, "TimeStart"))
		{
			if (fscanf(f, "%s", s2) == EOF) return false;
			entry.time_start = atof(s2);
		}
		else
		{
			fsetpos(f, &pos2);
		}
	}

	if (entry.ids.empty()) return false;

	userdef_elem_entries.push_back(entry);
	return true;
}

bool Monitor::ReadElemUserDefParams_B(const char* s, UserDefElemMonitorParams& p)
{
	// double (escalar)
	if      (!strcmp(s, "length"))           p.length           = true;
	else if (!strcmp(s, "jacobian"))         p.jacobian         = true;
	else if (!strcmp(s, "alpha1"))           p.alpha1           = true;
	else if (!strcmp(s, "p0i"))              p.p0i              = true;
	else if (!strcmp(s, "p0e"))              p.p0e              = true;
	else if (!strcmp(s, "rhoi"))             p.rhoi             = true;
	else if (!strcmp(s, "rhoe"))             p.rhoe             = true;
	else if (!strcmp(s, "Aint"))             p.Aint             = true;
	else if (!strcmp(s, "rho_adt"))          p.rho_adt          = true;
	else if (!strcmp(s, "rho_adn"))          p.rho_adn          = true;
	else if (!strcmp(s, "load_multiplier"))  p.load_multiplier  = true;
	else if (!strcmp(s, "alpha_escalar_delta")) p.alpha_escalar_delta = true;
	else if (!strcmp(s, "g"))                p.g                = true;
	else if (!strcmp(s, "signal_t"))         p.signal_t         = true;
	else if (!strcmp(s, "Cdt"))              p.Cdt              = true;
	else if (!strcmp(s, "Cdn"))              p.Cdn              = true;
	else if (!strcmp(s, "Aext"))             p.Aext             = true;
	else if (!strcmp(s, "rho_f"))            p.rho_f            = true;
	else if (!strcmp(s, "depth"))            p.depth            = true;
	else if (!strcmp(s, "Un_"))              p.Un_              = true;
	else if (!strcmp(s, "un_"))              p.un_scalar        = true;
	else if (!strcmp(s, "Ut_"))              p.Ut_              = true;
	else if (!strcmp(s, "ut_"))              p.ut_scalar        = true;
	else if (!strcmp(s, "C1t"))              p.C1t              = true;
	else if (!strcmp(s, "C1n"))              p.C1n              = true;
	else if (!strcmp(s, "l_factor"))         p.l_factor         = true;
	else if (!strcmp(s, "mult"))             p.mult             = true;
	else if (!strcmp(s, "t1"))               p.t1               = true;
	else if (!strcmp(s, "t2"))               p.t2               = true;
	// Matrix* (matriz unica)
	else if (!strcmp(s, "e3r"))              p.e3r              = true;
	else if (!strcmp(s, "e3rg"))             p.e3rg             = true;
	else if (!strcmp(s, "omega_rb"))         p.omega_rb         = true;
	else if (!strcmp(s, "r_inst"))           p.r_inst           = true;
	else if (!strcmp(s, "i_loading"))        p.i_loading        = true;
	else if (!strcmp(s, "e_loading"))        p.e_loading        = true;
	else if (!strcmp(s, "P_loading"))        p.P_loading        = true;
	else if (!strcmp(s, "inertial_loading")) p.inertial_loading = true;
	else if (!strcmp(s, "morison_loading"))  p.morison_loading  = true;
	else if (!strcmp(s, "transform3"))       p.transform3       = true;
	else if (!strcmp(s, "Mr"))               p.Mr_mat           = true;
	else if (!strcmp(s, "Jr"))               p.Jr_mat           = true;
	else if (!strcmp(s, "stiffness"))        p.stiffness_mat    = true;
	else if (!strcmp(s, "D"))                p.D_mat            = true;
	else if (!strcmp(s, "I3"))               p.I3               = true;
	else if (!strcmp(s, "B1"))               p.B1               = true;
	else if (!strcmp(s, "Qtransp"))          p.Qtransp          = true;
	else if (!strcmp(s, "B2"))               p.B2               = true;
	else if (!strcmp(s, "B2temp"))           p.B2temp           = true;
	else if (!strcmp(s, "constitutive_stiffness")) p.constitutive_stiffness = true;
	else if (!strcmp(s, "geometric_stiffness"))    p.geometric_stiffness    = true;
	else if (!strcmp(s, "loading_stiffness"))      p.loading_stiffness      = true;
	else if (!strcmp(s, "mass"))             p.mass             = true;
	else if (!strcmp(s, "damping"))          p.damping          = true;
	else if (!strcmp(s, "mass_modal"))       p.mass_modal       = true;
	else if (!strcmp(s, "damping_modal"))    p.damping_modal    = true;
	else if (!strcmp(s, "damping_loading"))  p.damping_loading  = true;
	else if (!strcmp(s, "rayleigh_damping")) p.rayleigh_damping = true;
	else if (!strcmp(s, "transform"))        p.transform        = true;
	else if (!strcmp(s, "V_alpha_dz_n"))     p.V_alpha_dz_n     = true;
	else if (!strcmp(s, "V_alpha_m"))        p.V_alpha_m        = true;
	else if (!strcmp(s, "d_V_dalpha_apha_m"))p.d_V_dalpha_apha_m= true;
	else if (!strcmp(s, "G_d_u_alpha"))      p.G_d_u_alpha      = true;
	else if (!strcmp(s, "G_d_u_alpha_transp"))p.G_d_u_alpha_transp= true;
	else if (!strcmp(s, "G_alpha_alpha"))    p.G_alpha_alpha    = true;
	else if (!strcmp(s, "G_alpha_d_alpha"))  p.G_alpha_d_alpha  = true;
	else if (!strcmp(s, "G_alpha_d_alpha_transp")) p.G_alpha_d_alpha_transp = true;
	else if (!strcmp(s, "B"))                p.B_mat            = true;
	else if (!strcmp(s, "G"))                p.G_mat            = true;
	else if (!strcmp(s, "K1ua"))             p.K1ua             = true;
	else if (!strcmp(s, "K1aa"))             p.K1aa             = true;
	else if (!strcmp(s, "K2ua"))             p.K2ua             = true;
	else if (!strcmp(s, "K2au"))             p.K2au             = true;
	else if (!strcmp(s, "Kext"))             p.Kext             = true;
	else if (!strcmp(s, "O1"))               p.O1               = true;
	else if (!strcmp(s, "Lambda"))           p.Lambda           = true;
	// vector<double>
	else if (!strcmp(s, "external_pressure"))p.external_pressure= true;
	// LengthPoints
	else if (!strcmp(s, "length_points"))    p.length_points    = true;
	else return false;
	return true;
}

bool Monitor::ReadElemUserDefParams(FILE *f, UserDefElemMonitorParams& p)
{
	char s[1000];
	fpos_t pos;
	bool found = false;

	while (true)
	{
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF) break;

		bool is_param = true;
		// Matrix** (2 GPs)
		if      (!strcmp(s, "sigma_r"))          p.sigma_r          = true;
		// Componentes individuais de sigma_r
		else if (!strcmp(s, "nx1r"))             p.nx1r             = true;
		else if (!strcmp(s, "ny1r"))             p.ny1r             = true;
		else if (!strcmp(s, "nz1r"))             p.nz1r             = true;
		else if (!strcmp(s, "mx1r"))             p.mx1r             = true;
		else if (!strcmp(s, "my1r"))             p.my1r             = true;
		else if (!strcmp(s, "mz1r"))             p.mz1r             = true;
		else if (!strcmp(s, "nx2r"))             p.nx2r             = true;
		else if (!strcmp(s, "ny2r"))             p.ny2r             = true;
		else if (!strcmp(s, "nz2r"))             p.nz2r             = true;
		else if (!strcmp(s, "mx2r"))             p.mx2r             = true;
		else if (!strcmp(s, "my2r"))             p.my2r             = true;
		else if (!strcmp(s, "mz2r"))             p.mz2r             = true;
		else if (!strcmp(s, "eta_r"))            p.eta_r            = true;
		else if (!strcmp(s, "kappa_r"))          p.kappa_r          = true;
		else if (!strcmp(s, "epsilon_r"))        p.epsilon_r        = true;
		else if (!strcmp(s, "n_r"))              p.n_r              = true;
		else if (!strcmp(s, "m_r"))              p.m_r              = true;
		else if (!strcmp(s, "n_global"))         p.n_global         = true;
		else if (!strcmp(s, "m_global"))         p.m_global         = true;
		else if (!strcmp(s, "alpha_delta"))      p.alpha_delta      = true;
		else if (!strcmp(s, "d_alpha_delta"))    p.d_alpha_delta    = true;
		else if (!strcmp(s, "u_delta"))          p.u_delta          = true;
		else if (!strcmp(s, "d_u_delta"))        p.d_u_delta        = true;
		else if (!strcmp(s, "A_delta"))          p.A_delta          = true;
		else if (!strcmp(s, "Q_delta"))          p.Q_delta          = true;
		else if (!strcmp(s, "d_z"))              p.d_z              = true;
		else if (!strcmp(s, "d_Z"))              p.d_Z              = true;
		else if (!strcmp(s, "Xi_delta"))         p.Xi_delta         = true;
		else if (!strcmp(s, "d_A_delta"))        p.d_A_delta        = true;
		else if (!strcmp(s, "d_Xi_delta"))       p.d_Xi_delta       = true;
		// Fun��es de forma (Matrix** nos GPs)
		else if (!strcmp(s, "N_gp"))             p.N_gp             = true;
		else if (!strcmp(s, "deltaN"))           p.deltaN_gp        = true;
		// Correnteza mar�tima (Matrix** nos GPs)
		else if (!strcmp(s, "e3ip"))             p.e3ip             = true;
		else if (!strcmp(s, "zi"))               p.zi               = true;
		else if (!strcmp(s, "vel"))              p.vel              = true;
		else if (!strcmp(s, "velr"))             p.velr             = true;
		else if (!strcmp(s, "element_vel"))      p.element_vel      = true;
		else if (!strcmp(s, "ut"))               p.ut_gp            = true;
		else if (!strcmp(s, "un"))               p.un_gp            = true;
		else if (!strcmp(s, "d_e3_d_alpha"))     p.d_e3_d_alpha     = true;
		else if (!strcmp(s, "Lt"))               p.Lt               = true;
		else if (!strcmp(s, "Ln"))               p.Ln               = true;
		else if (!strcmp(s, "L_u_alpha"))        p.L_u_alpha        = true;
		else if (!strcmp(s, "f_current"))        p.f_current        = true;
		else if (!strcmp(s, "L"))                p.L_gp             = true;
		else if (!strcmp(s, "t_e"))              p.t_e              = true;
		else if (!strcmp(s, "n_e"))              p.n_e              = true;
		else if (!strcmp(s, "vtr"))              p.vtr              = true;
		else if (!strcmp(s, "vnr"))              p.vnr              = true;
		else if (!strcmp(s, "Mdt"))              p.Mdt              = true;
		else if (!strcmp(s, "Mdn"))              p.Mdn              = true;
		else if (!strcmp(s, "Md2"))              p.Md2              = true;
		// Pipe load (Matrix** nos GPs)
		else if (!strcmp(s, "kip"))              p.kip              = true;
		else if (!strcmp(s, "temp_f"))           p.temp_f           = true;
		else if (!strcmp(s, "temp_m"))           p.temp_m           = true;
		else if (!strcmp(s, "temp_l"))           p.temp_l           = true;
		else if (!strcmp(s, "Kip"))              p.Kip              = true;
		else if (!strcmp(s, "E3ip"))             p.E3ip             = true;
		else if (!strcmp(s, "UpsilonN"))         p.UpsilonN         = true;
		// Din�mica (Matrix** nos GPs)
		else if (!strcmp(s, "Xi_dot"))           p.Xi_dot           = true;
		else if (!strcmp(s, "Mip"))              p.Mip_gp           = true;
		else if (!strcmp(s, "Jip"))              p.Jip_gp           = true;
		else if (!strcmp(s, "M"))                p.M_dyn            = true;
		else if (!strcmp(s, "Md1"))              p.Md1              = true;
		else if (!strcmp(s, "omega_ip"))         p.omega_ip         = true;
		else if (!strcmp(s, "temp_if_l"))        p.temp_if_l        = true;
		else if (!strcmp(s, "f_pressure_1"))     p.f_pressure_1     = true;
		else if (!strcmp(s, "f_pressure_2"))     p.f_pressure_2     = true;
		else if (!strcmp(s, "alpha_dot"))        p.alpha_dot        = true;
		// double* (2 GPs)
		else if (!strcmp(s, "pressure"))         p.pressure         = true;
		else if (!strcmp(s, "rho"))              p.rho              = true;
		else if (!strcmp(s, "temperature"))      p.temperature      = true;
		else if (!strcmp(s, "flow_velocity"))    p.flow_velocity    = true;
		else if (!strcmp(s, "flow_rate"))        p.flow_rate        = true;
		else if (!strcmp(s, "internal_pressure"))p.internal_pressure= true;
		else if (!strcmp(s, "N1"))               p.N1               = true;
		else if (!strcmp(s, "N2"))               p.N2               = true;
		else if (!strcmp(s, "N3"))               p.N3               = true;
		else if (!strcmp(s, "dN1"))              p.dN1              = true;
		else if (!strcmp(s, "dN2"))              p.dN2              = true;
		else if (!strcmp(s, "dN3"))              p.dN3              = true;
		else if (!strcmp(s, "csi"))              p.csi              = true;
		else if (!ReadElemUserDefParams_B(s, p)) is_param = false;

		if (!is_param)
		{
			fsetpos(f, &pos);
			break;
		}
		found = true;
	}
	return found;
}
void Monitor::Write(FILE *f)
{
	fprintf(f, "Monitor\tSample\t%d\n",sample);
	if (monitor_nodes_exists)
	{
		fprintf(f, "MonitorNodes\t");
		for (int i = 0; i < nodes.size(); i++)
			fprintf(f, "%d\t", nodes[i]);
		fprintf(f, "\n");
	}
	if (monitor_elements_exists)
	{
		fprintf(f, "MonitorElements\t");
		for (int i = 0; i < elements.size(); i++)
			fprintf(f, "%d\t", elements[i]);
		fprintf(f, "\n");
	}
	if (monitor_contacts_exists)
	{
		fprintf(f, "MonitorContacts\t");
		for (int i = 0; i < contacts.size(); i++)
			fprintf(f, "%d\t", contacts[i]);
		fprintf(f, "\n");
	}
	if (monitor_node_sets_exists)
	{
		fprintf(f, "MonitorNodeSets\t");
		for (int i = 0; i < node_sets.size(); i++)
			fprintf(f, "%d\t", node_sets[i]);
		fprintf(f, "\n");
	}
	if (monitor_particles_exists)
	{
		fprintf(f, "MonitorParticles\t");
		for (int i = 0; i < particles.size(); i++)
			fprintf(f, "%d\t", particles[i]);
		fprintf(f, "\n");
	}
	if (monitor_global_exists)
	{
		fprintf(f, "MonitorGlobal\t");
		fprintf(f, "\n");
	}

	if (monitor_special_constraints_exists)
	{
		fprintf(f, "MonitorSpecialConstraints\t");
		for (int i = 0; i < special_constraints.size(); i++)
			fprintf(f, "%d\t", special_constraints[i]);
		fprintf(f, "\n");
	}

	for (int e = 0; e < (int)userdef_entries.size(); e++)
	{
		const UserDefMonitorEntry& entry = userdef_entries[e];
		fprintf(f, "MonitorNodesUserDefined\t");
		fprintf(f, entry.use_node_set ? "NodeSet\t" : "Nodes\t");
		for (int id : entry.ids) fprintf(f, "%d\t", id);
		fprintf(f, "Parameters\t");
		const UserDefMonitorParams& p = entry.params;
		if (p.X)   fprintf(f, "X\t");
		if (p.Y)   fprintf(f, "Y\t");
		if (p.Z)   fprintf(f, "Z\t");
		if (p.dX)  fprintf(f, "dX\t");
		if (p.dY)  fprintf(f, "dY\t");
		if (p.dZ)  fprintf(f, "dZ\t");
		if (p.ddX) fprintf(f, "ddX\t");
		if (p.ddY) fprintf(f, "ddY\t");
		if (p.ddZ) fprintf(f, "ddZ\t");
		if (p.FX)  fprintf(f, "FX\t");
		if (p.FY)  fprintf(f, "FY\t");
		if (p.FZ)  fprintf(f, "FZ\t");
		if (p.MX)  fprintf(f, "MX\t");
		if (p.MY)  fprintf(f, "MY\t");
		if (p.MZ)  fprintf(f, "MZ\t");
		if (entry.time_start != 0.0) fprintf(f, "TimeStart\t%.6e\t", entry.time_start);
		fprintf(f, "\n");
	}

	for (int e = 0; e < (int)userdef_elem_entries.size(); e++)
	{
		const UserDefElemMonitorEntry& entry = userdef_elem_entries[e];
		fprintf(f, "MonitorElementsUserDefined\t");
		fprintf(f, entry.use_element_set ? "ElementSet\t" : "Elements\t");
		for (int id : entry.ids) fprintf(f, "%d\t", id);
		fprintf(f, "Parameters\t");
		const UserDefElemMonitorParams& p = entry.params;
		// Matrix** (2 GPs)
		if (p.sigma_r)          fprintf(f, "sigma_r\t");
		// Componentes individuais de sigma_r
		if (p.nx1r) fprintf(f, "nx1r\t"); if (p.ny1r) fprintf(f, "ny1r\t"); if (p.nz1r) fprintf(f, "nz1r\t");
		if (p.mx1r) fprintf(f, "mx1r\t"); if (p.my1r) fprintf(f, "my1r\t"); if (p.mz1r) fprintf(f, "mz1r\t");
		if (p.nx2r) fprintf(f, "nx2r\t"); if (p.ny2r) fprintf(f, "ny2r\t"); if (p.nz2r) fprintf(f, "nz2r\t");
		if (p.mx2r) fprintf(f, "mx2r\t"); if (p.my2r) fprintf(f, "my2r\t"); if (p.mz2r) fprintf(f, "mz2r\t");
		if (p.eta_r)            fprintf(f, "eta_r\t");
		if (p.kappa_r)          fprintf(f, "kappa_r\t");
		if (p.epsilon_r)        fprintf(f, "epsilon_r\t");
		if (p.n_r)              fprintf(f, "n_r\t");
		if (p.m_r)              fprintf(f, "m_r\t");
		if (p.n_global)         fprintf(f, "n_global\t");
		if (p.m_global)         fprintf(f, "m_global\t");
		if (p.alpha_delta)      fprintf(f, "alpha_delta\t");
		if (p.d_alpha_delta)    fprintf(f, "d_alpha_delta\t");
		if (p.u_delta)          fprintf(f, "u_delta\t");
		if (p.d_u_delta)        fprintf(f, "d_u_delta\t");
		if (p.A_delta)          fprintf(f, "A_delta\t");
		if (p.Q_delta)          fprintf(f, "Q_delta\t");
		if (p.d_z)              fprintf(f, "d_z\t");
		if (p.omega_ip)         fprintf(f, "omega_ip\t");
		if (p.temp_if_l)        fprintf(f, "temp_if_l\t");
		if (p.f_pressure_1)     fprintf(f, "f_pressure_1\t");
		if (p.f_pressure_2)     fprintf(f, "f_pressure_2\t");
		if (p.alpha_dot)        fprintf(f, "alpha_dot\t");
		// double* (2 GPs)
		if (p.pressure)         fprintf(f, "pressure\t");
		if (p.rho)              fprintf(f, "rho\t");
		if (p.temperature)      fprintf(f, "temperature\t");
		if (p.flow_velocity)    fprintf(f, "flow_velocity\t");
		if (p.flow_rate)        fprintf(f, "flow_rate\t");
		if (p.internal_pressure)fprintf(f, "internal_pressure\t");
		// double (escalar)
		if (p.length)           fprintf(f, "length\t");
		if (p.jacobian)         fprintf(f, "jacobian\t");
		if (p.alpha1)           fprintf(f, "alpha1\t");
		if (p.p0i)              fprintf(f, "p0i\t");
		if (p.p0e)              fprintf(f, "p0e\t");
		if (p.rhoi)             fprintf(f, "rhoi\t");
		if (p.rhoe)             fprintf(f, "rhoe\t");
		if (p.Aint)             fprintf(f, "Aint\t");
		if (p.rho_adt)          fprintf(f, "rho_adt\t");
		if (p.rho_adn)          fprintf(f, "rho_adn\t");
		if (p.load_multiplier)  fprintf(f, "load_multiplier\t");
		// Matrix* (matriz unica)
		if (p.e3r)              fprintf(f, "e3r\t");
		if (p.e3rg)             fprintf(f, "e3rg\t");
		if (p.omega_rb)         fprintf(f, "omega_rb\t");
		if (p.r_inst)           fprintf(f, "r_inst\t");
		if (p.i_loading)        fprintf(f, "i_loading\t");
		if (p.e_loading)        fprintf(f, "e_loading\t");
		if (p.P_loading)        fprintf(f, "P_loading\t");
		if (p.inertial_loading) fprintf(f, "inertial_loading\t");
		if (p.morison_loading)  fprintf(f, "morison_loading\t");
		if (p.transform3)       fprintf(f, "transform3\t");
		if (p.Mr_mat)           fprintf(f, "Mr\t");
		if (p.Jr_mat)           fprintf(f, "Jr\t");
		if (p.stiffness_mat)    fprintf(f, "stiffness\t");
		if (p.D_mat)            fprintf(f, "D\t");
		// vector<double>
		if (p.external_pressure)fprintf(f, "external_pressure\t");
		// LengthPoints
		if (p.length_points)    fprintf(f, "length_points\t");
		if (entry.time_start != 0.0) fprintf(f, "TimeStart\t%.6e\t", entry.time_start);
		fprintf(f, "\n");
	}	
}

void Monitor::StartMonitor()
{
	//Abre os arquivos dos n�s
	char name[200];
	for (int i = 0; i < nodes.size(); i++)
	{
		int node = nodes[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(node, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_node_");
		strcat(name, number);
		strcat(name, ".txt");
		f_nodes[i] = fopen(name, "w");
		//fprintf(f_nodes[i], "NODE\t%d\n",node);
	}
	
	//Abre os arquivos dos elementos
	for (int i = 0; i < elements.size(); i++)
	{
		int element = elements[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(element, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_element_");
		strcat(name, number);
		strcat(name, ".txt");
		f_elements[i] = fopen(name, "w");
		//fprintf(f_elements[i], "ELEMENT\t%d\t%s\n", element,db.elements[element-1]->type_name);
	}

	//Abre os arquivos dos contatos
	for (int i = 0; i < contacts.size(); i++)
	{
		int contact = contacts[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(contact, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_contact_");
		strcat(name, number);
		strcat(name, ".txt");
		f_contacts[i] = fopen(name, "w");
		//fprintf(f_contacts[i], "CONTACT\t%d\t%s\n", contact, db.contacts[contact-1]->type_name);
	}

	//Abre os arquivos dos node sets
	for (int i = 0; i < node_sets.size(); i++)
	{
		int set = node_sets[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(set, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_nodeset_");
		strcat(name, number);
		strcat(name, ".txt");
		f_node_sets[i] = fopen(name, "w");
		//fprintf(f_node_sets[i], "NODESET\t%d\n", set);
	}

	//Abre os arquivos das particulas
	for (int i = 0; i < particles.size(); i++)
	{
		int particle = particles[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(particle, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_particle_");
		strcat(name, number);
		strcat(name, ".txt");
		f_particles[i] = fopen(name, "w");
		//fprintf(f_elements[i], "PARTICLE\t%d\t%s\n", element,db.particles[particle-1]->type_name);
	}

	//Abre o arquivo do monitor global
	strcpy(name, db.folder_name);
	strcat(name, "monitors/");
	_mkdir(name);
	strcat(name, "monitor_global");
	strcat(name, ".txt");
	f_global[0] = fopen(name, "w");

	//Abre os arquivos dos special constraints
	for (int i = 0; i < special_constraints.size(); i++)
	{
		int sc = special_constraints[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(sc, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_special_constraint_");
		strcat(name, number);
		strcat(name, ".txt");
		f_special_constraints[i] = fopen(name, "w");
		//fprintf(f_elements[i], "SPECIALCONSTRAINT\t%d\t%s\n", sc,db.special_constraints[sc-1]->type_name);
	}
	
// Abre os arquivos dos monitores de n�s definidos pelo usu�rio
	for (int e = 0; e < (int)userdef_entries.size(); e++)
	{
		UserDefMonitorEntry& entry = userdef_entries[e];

		// Resolver IDs para node IDs concretos
		entry.resolved_node_ids.clear();
		if (entry.use_node_set)
		{
			for (int set_id : entry.ids)
			{
				NodeSet* ns = db.node_sets[set_id - 1];
				for (int k = 0; k < ns->n_nodes; k++)
					entry.resolved_node_ids.push_back(ns->node_list[k]);
			}
		}
		else
		{
			entry.resolved_node_ids = entry.ids;
		}
		// Inicializar node_params, node_time_starts e node_first_record para cada n�
		int n = (int)entry.resolved_node_ids.size();
		entry.node_params.assign(n, entry.params);
		entry.node_time_starts.assign(n, entry.time_start);
		entry.node_first_record.assign(n, true);
	}

	// Consolidar n�s duplicados entre entradas: OR dos params por n�
	for (int e1 = 0; e1 < (int)userdef_entries.size(); e1++)
	{
		for (int e2 = e1 + 1; e2 < (int)userdef_entries.size(); )
		{
			UserDefMonitorEntry& en1 = userdef_entries[e1];
			UserDefMonitorEntry& en2 = userdef_entries[e2];

			vector<int>              remaining_ids;
			vector<UserDefMonitorParams> remaining_params;
			vector<double>           remaining_time_starts;
			vector<bool>             remaining_first_record;
			for (int j = 0; j < (int)en2.resolved_node_ids.size(); j++)
			{
				int node_id = en2.resolved_node_ids[j];
				auto it = find(en1.resolved_node_ids.begin(),
				               en1.resolved_node_ids.end(), node_id);
				if (it != en1.resolved_node_ids.end())
				{
					int idx1 = (int)(it - en1.resolved_node_ids.begin());
					UserDefMonitorParams&       p1 = en1.node_params[idx1];
					const UserDefMonitorParams& p2 = en2.node_params[j];
					p1.X   |= p2.X;   p1.Y   |= p2.Y;   p1.Z   |= p2.Z;
					p1.dX  |= p2.dX;  p1.dY  |= p2.dY;  p1.dZ  |= p2.dZ;
					p1.ddX |= p2.ddX; p1.ddY |= p2.ddY; p1.ddZ |= p2.ddZ;
					p1.FX  |= p2.FX;  p1.FY  |= p2.FY;  p1.FZ  |= p2.FZ;
					p1.MX  |= p2.MX;  p1.MY  |= p2.MY;  p1.MZ  |= p2.MZ;
					// time_start: usar o menor entre as duas entradas
					if (en2.node_time_starts[j] < en1.node_time_starts[idx1])
						en1.node_time_starts[idx1] = en2.node_time_starts[j];
				}
				else
				{
					remaining_ids.push_back(node_id);
					remaining_params.push_back(en2.node_params[j]);
					remaining_time_starts.push_back(en2.node_time_starts[j]);
					remaining_first_record.push_back(true);
				}
			}
			en2.resolved_node_ids = remaining_ids;
			en2.node_params       = remaining_params;
			en2.node_time_starts  = remaining_time_starts;
			en2.node_first_record = remaining_first_record;

			if (en2.resolved_node_ids.empty())
				userdef_entries.erase(userdef_entries.begin() + e2);
			else
				e2++;
		}
	}

	// === MonitorElementsUserDefined: resolver IDs, consolidar duplicatas, abrir arquivos ===
	for (int e = 0; e < (int)userdef_elem_entries.size(); e++)
	{
		UserDefElemMonitorEntry& entry = userdef_elem_entries[e];

		entry.resolved_elem_ids.clear();
		if (entry.use_element_set)
		{
			for (int set_id : entry.ids)
			{
				ElementSet* es = db.element_sets[set_id - 1];
				for (int k = 0; k < es->n_el; k++)
					entry.resolved_elem_ids.push_back(es->el_list[k]);
			}
		}
		else
		{
			entry.resolved_elem_ids = entry.ids;
		}

		int n = (int)entry.resolved_elem_ids.size();
		entry.elem_params.assign(n, entry.params);
		entry.elem_time_starts.assign(n, entry.time_start);
		entry.elem_first_record.assign(n, true);
	}

	// Consolidar elementos duplicados entre entradas: OR dos params, min do time_start
	for (int e1 = 0; e1 < (int)userdef_elem_entries.size(); e1++)
	{
		for (int e2 = e1 + 1; e2 < (int)userdef_elem_entries.size(); )
		{
			UserDefElemMonitorEntry& en1 = userdef_elem_entries[e1];
			UserDefElemMonitorEntry& en2 = userdef_elem_entries[e2];

			vector<int>                  remaining_ids;
			vector<UserDefElemMonitorParams> remaining_params;
			vector<double>               remaining_time_starts;
			vector<bool>                 remaining_first_record;

			for (int j = 0; j < (int)en2.resolved_elem_ids.size(); j++)
			{
				int elem_id = en2.resolved_elem_ids[j];
				auto it = find(en1.resolved_elem_ids.begin(), en1.resolved_elem_ids.end(), elem_id);
				if (it != en1.resolved_elem_ids.end())
				{
					int idx1 = (int)(it - en1.resolved_elem_ids.begin());
					UserDefElemMonitorParams&       p1 = en1.elem_params[idx1];
					const UserDefElemMonitorParams& p2 = en2.elem_params[j];
					p1.sigma_r       |= p2.sigma_r;
					p1.nx1r |= p2.nx1r; p1.ny1r |= p2.ny1r; p1.nz1r |= p2.nz1r;
					p1.mx1r |= p2.mx1r; p1.my1r |= p2.my1r; p1.mz1r |= p2.mz1r;
					p1.nx2r |= p2.nx2r; p1.ny2r |= p2.ny2r; p1.nz2r |= p2.nz2r;
					p1.mx2r |= p2.mx2r; p1.my2r |= p2.my2r; p1.mz2r |= p2.mz2r;
					p1.eta_r         |= p2.eta_r;
					p1.kappa_r       |= p2.kappa_r;
					p1.epsilon_r     |= p2.epsilon_r;
					p1.n_r           |= p2.n_r;
					p1.m_r           |= p2.m_r;
					p1.n_global      |= p2.n_global;
					p1.m_global      |= p2.m_global;
					p1.alpha_delta   |= p2.alpha_delta;
					p1.d_alpha_delta |= p2.d_alpha_delta;
					p1.u_delta       |= p2.u_delta;
					p1.d_u_delta     |= p2.d_u_delta;
					p1.pressure      |= p2.pressure;
					p1.rho           |= p2.rho;
					p1.temperature   |= p2.temperature;
					p1.flow_velocity |= p2.flow_velocity;
					p1.flow_rate     |= p2.flow_rate;
					p1.length        |= p2.length;
					if (en2.elem_time_starts[j] < en1.elem_time_starts[idx1])
						en1.elem_time_starts[idx1] = en2.elem_time_starts[j];
				}
				else
				{
					remaining_ids.push_back(elem_id);
					remaining_params.push_back(en2.elem_params[j]);
					remaining_time_starts.push_back(en2.elem_time_starts[j]);
					remaining_first_record.push_back(true);
				}
			}
			en2.resolved_elem_ids  = remaining_ids;
			en2.elem_params        = remaining_params;
			en2.elem_time_starts   = remaining_time_starts;
			en2.elem_first_record  = remaining_first_record;

			if (en2.resolved_elem_ids.empty())
				userdef_elem_entries.erase(userdef_elem_entries.begin() + e2);
			else
				e2++;
		}
	}

	// Abrir um arquivo por elemento resolvido
	for (int e = 0; e < (int)userdef_elem_entries.size(); e++)
	{
		UserDefElemMonitorEntry& entry = userdef_elem_entries[e];
		entry.files.resize(entry.resolved_elem_ids.size(), nullptr);
		for (int k = 0; k < (int)entry.resolved_elem_ids.size(); k++)
		{
			int elem_id = entry.resolved_elem_ids[k];
			char num_e[20];
			strcpy(name, db.folder_name);
			strcat(name, "monitors/");
			_mkdir(name);
			_itoa(elem_id, num_e, 10);
			strcat(name, "monitor_userdef_element_");
			strcat(name, num_e);
			strcat(name, ".txt");
			entry.files[k] = fopen(name, "w");
		}
	}

	// Abrir um arquivo por n� resolvido
	for (int e = 0; e < (int)userdef_entries.size(); e++)
	{
		UserDefMonitorEntry& entry = userdef_entries[e];
		entry.files.resize(entry.resolved_node_ids.size(), nullptr);
		for (int k = 0; k < (int)entry.resolved_node_ids.size(); k++)
		{
			int node_id = entry.resolved_node_ids[k];
			char num_n[20];
			strcpy(name, db.folder_name);
			strcat(name, "monitors/");
			_mkdir(name);
			_itoa(node_id, num_n, 10);
			strcat(name, "monitor_userdef_node_");
			strcat(name, num_n);
			strcat(name, ".txt");
			entry.files[k] = fopen(name, "w");
		}
	}
}

void Monitor::UpdateMonitor(double time)
{
	//Atualiza��es a serem realizadas - salvas nos arquivos do monitor
	for (int i = 0; i < nodes.size(); i++)
		db.nodes[nodes[i] - 1]->WriteMonitor(f_nodes[i], first_record,time);
	for (int i = 0; i < elements.size(); i++)
		db.elements[elements[i] - 1]->WriteMonitor(f_elements[i], first_record, time);
	for (int i = 0; i < contacts.size(); i++)
		db.contacts[contacts[i] - 1]->WriteMonitor(f_contacts[i], first_record, time);
	for (int i = 0; i < node_sets.size(); i++)
		db.node_sets[node_sets[i] - 1]->WriteMonitor(f_node_sets[i], first_record, time);
	for (int i = 0; i < particles.size(); i++)
		db.particles[particles[i] - 1]->WriteMonitor(f_particles[i], first_record, time);
	UpdateGlobalMonitor(time);
	for (int i = 0; i < special_constraints.size(); i++)
		db.special_constraints[special_constraints[i] - 1]->WriteMonitor(f_special_constraints[i], first_record, time);
	// Monitores de n�s definidos pelo usu�rio
	for (auto& entry : userdef_entries)
		for (int k = 0; k < (int)entry.resolved_node_ids.size(); k++)
		{
			if (time < entry.node_time_starts[k]) continue;
			db.nodes[entry.resolved_node_ids[k] - 1]->WriteMonitorUserDef(
				entry.files[k], entry.node_first_record[k], time, entry.node_params[k]);
			entry.node_first_record[k] = false;
		}

	// Monitores de elementos definidos pelo usu�rio
	for (auto& entry : userdef_elem_entries)
		for (int k = 0; k < (int)entry.resolved_elem_ids.size(); k++)
		{
			if (time < entry.elem_time_starts[k]) continue;
			static_cast<Pipe_1*>(db.elements[entry.resolved_elem_ids[k] - 1])->WriteMonitorUserDef(
				entry.files[k], entry.elem_first_record[k], time, entry.elem_params[k]);
			entry.elem_first_record[k] = false;
		}
	first_record = false;
}

void Monitor::EndMonitor()
{
	for (int i = 0; i < nodes.size(); i++)
		fclose(f_nodes[i]);
	for (int i = 0; i < elements.size(); i++)
		fclose(f_elements[i]);
	for (int i = 0; i < contacts.size(); i++)
		fclose(f_contacts[i]);
	for (int i = 0; i < node_sets.size(); i++)
		fclose(f_node_sets[i]);
	for (int i = 0; i < special_constraints.size(); i++)
		fclose(f_special_constraints[i]);
	for (int i = 0; i < particles.size(); i++)
		fclose(f_particles[i]);
	fclose(f_global[0]);
		for (auto& entry : userdef_entries)
	{
		for (FILE* fp : entry.files)
			if (fp) fclose(fp);
		entry.files.clear();
	}

	for (auto& entry : userdef_elem_entries)
	{
		for (FILE* fp : entry.files)
			if (fp) fclose(fp);
		entry.files.clear();
	}
}

void Monitor::UpdateGlobalMonitor(double time)
{
	//Cabe�alho
	if (first_record == true)
	{
		fprintf(f_global[0], "TIME\t");
		fprintf(f_global[0], "KIN_EN\t");
		fprintf(f_global[0], "STRAIN_EN\t");
		fprintf(f_global[0], "GRAV_EN\t");
		fprintf(f_global[0], "MECH_EN\t");
		/*fprintf(f_global[0], "LIN_MOM_X\t");
		fprintf(f_global[0], "LIN_MOM_Y\t");
		fprintf(f_global[0], "LIN_MOM_Z\t");
		fprintf(f_global[0], "ANG_MOM_X\t");
		fprintf(f_global[0], "ANG_MOM_Y\t");
		fprintf(f_global[0], "ANG_MOM_Z\t");*/
		if (db.gcs_exist)
		{
			fprintf(f_global[0], "CONTACT_ENERGY\t");
			fprintf(f_global[0], "MONITORED_PP\t");
			fprintf(f_global[0], "ACTIVE_PP\t");
			fprintf(f_global[0], "MONITORED_PB\t");
			fprintf(f_global[0], "ACTIVE_PB\t");
			fprintf(f_global[0], "MONITORED_BOBO\t");
			fprintf(f_global[0], "ACTIVE_BOBO\t");
			fprintf(f_global[0], "MONITORED_PBO\t");
			fprintf(f_global[0], "ACTIVE_PBO\t");

			if (print_times)
			{
				fprintf(f_global[0], "DURATION_VERLET\t");
				fprintf(f_global[0], "DURATION_LCELLS\t");
				fprintf(f_global[0], "DURATION_COLDETECTION\t");
				fprintf(f_global[0], "NUMBER_COLDETECTION\t");
				fprintf(f_global[0], "DURATION_MOUNTCONTACTS\t");
			}
		}
		
		fprintf(f_global[0], "TIMESTEP\t");
		fprintf(f_global[0], "F_CRIT\t");
		fprintf(f_global[0], "M_CRIT\t");
		fprintf(f_global[0], "DISP_CRIT\t");
		fprintf(f_global[0], "ROT_CRIT\t");
		fprintf(f_global[0], "\n");
	}
	//Informa��es a serem salvas
	fprintf(f_global[0], "%.6e\t",time);

	//Energia cinetica, momento linear e momento angular das particulas
	double kin = 0;
	double strain = 0;
	double grav = 0;
	/*Matrix linear_momentum(3);
	Matrix angular_momentum(3);*/
	if (db.particles_exist)
	{
		for (int i = 0; i < db.number_particles; i++)
		{
			kin += db.particles[i]->kinetic_energy;
			strain += db.particles[i]->strain_energy;
			grav += db.particles[i]->potential_g_energy;
			/*linear_momentum(0, 0) += db.particles[i]->linear_momentum[0];
			linear_momentum(1, 0) += db.particles[i]->linear_momentum[1];
			linear_momentum(2, 0) += db.particles[i]->linear_momentum[2];
			angular_momentum(0, 0) += db.particles[i]->angular_momentum_origin[0];
			angular_momentum(1, 0) += db.particles[i]->angular_momentum_origin[1];
			angular_momentum(2, 0) += db.particles[i]->angular_momentum_origin[2];*/
		}
			
	}
	if (db.elements_exist)
	{
		for (int i = 0; i < db.number_elements; i++)
		{
			kin += db.elements[i]->kinetic_energy;
			strain += db.elements[i]->strain_energy;
			grav += db.elements[i]->potential_gravitational_energy;
		}
	}
	fprintf(f_global[0], "%.6e\t", kin);
	fprintf(f_global[0], "%.6e\t", strain);
	fprintf(f_global[0], "%.6e\t", grav);
	fprintf(f_global[0], "%.6e\t", kin + strain + grav);
	/*fprintf(f_global[0], "%.6e\t%.6e\t%.6e\t", linear_momentum(0, 0), linear_momentum(1, 0), linear_momentum(2, 0));
	fprintf(f_global[0], "%.6e\t%.6e\t%.6e\t", angular_momentum(0, 0), angular_momentum(1, 0), angular_momentum(2, 0));*/

	if (db.gcs_exist)
	{
		fprintf(f_global[0], "%.6e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", db.gcs->EvaluateContactEnergy(), db.gcs->n_monitoring_PP, db.gcs->n_active_PP, db.gcs->n_monitoring_PB, db.gcs->n_active_PB,
			db.gcs->n_monitoring_BOBO, db.gcs->n_active_BOBO, db.gcs->n_monitoring_PBO, db.gcs->n_active_PBO);
		if (print_times)
		{
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_verlet);
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_linkedcells);
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_collision_detection);
			fprintf(f_global[0], "%d\t", (int)db.gcs->global_n_collisiondetection);
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_mount_contact);
		}
	}
		
	fprintf(f_global[0], "%.6e\t", db.current_time_step);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->force_criterion);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->moment_criterion);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->disp_criterion);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->rot_criterion);

	

	//Quebra de linha final
	fprintf(f_global[0], "\n");

	db.myprintf("\nSystem kinetic energy:\t%6e\n", kin);
	db.myprintf("System mechanical energy:\t%6e\n", kin + strain + grav);
}