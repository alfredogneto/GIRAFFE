#pragma once
#include <vector>
#include <stdio.h>
using namespace std;

// Par�metros selecion�veis por entrada do usu�rio para MonitorElementsUserDefined
// Tipos suportados:
//   double*    -> array de doubles nos 2 pontos de Gauss         -> 2 colunas por vari�vel
//   Matrix**   -> array de Matrix nos 2 pontos de Gauss          -> (m_lines*m_columns*2) colunas
//   double     -> escalar                                         -> 1 coluna
//   Matrix*    -> matriz �nica (qualquer dimens�o)                -> m_lines*m_columns colunas
//   vector<double> -> vetor de tamanho vari�vel                  -> N colunas (N = tamanho atual)
struct UserDefElemMonitorParams {

	// --- Matrix** (2 pontos de Gauss) ---
	// Componentes individuais de sigma_r (mesmos nomes do WriteMonitor padr�o)
	bool nx1r = false;	// sigma_r GP1 componente 0 (for�a normal x)
	bool ny1r = false;	// sigma_r GP1 componente 1 (for�a normal y)
	bool nz1r = false;	// sigma_r GP1 componente 2 (for�a normal z)
	bool mx1r = false;	// sigma_r GP1 componente 3 (momento x)
	bool my1r = false;	// sigma_r GP1 componente 4 (momento y)
	bool mz1r = false;	// sigma_r GP1 componente 5 (momento z)
	bool nx2r = false;	// sigma_r GP2 componente 0
	bool ny2r = false;	// sigma_r GP2 componente 1
	bool nz2r = false;	// sigma_r GP2 componente 2
	bool mx2r = false;	// sigma_r GP2 componente 3
	bool my2r = false;	// sigma_r GP2 componente 4
	bool mz2r = false;	// sigma_r GP2 componente 5
	// Deforma��es - sigma_r completo ou por componente
	bool sigma_r = false;	// tens�o (6x1 por GP)                -> 12 cols
	bool eta_r         = false;	// deforma��o extensional (3x1 por GP) ->  6 cols
	bool kappa_r       = false;	// curvatura (3x1 por GP)              ->  6 cols
	bool epsilon_r     = false;	// deforma��o total (6x1 por GP)       -> 12 cols
	// For�as e momentos internos (referencial local)
	bool n_r           = false;	// for�a interna local (3x1 por GP)    ->  6 cols
	bool m_r           = false;	// momento interno local (3x1 por GP)  ->  6 cols
	// For�as e momentos internos (referencial global)
	bool n_global      = false;	// for�a interna global "n" (3x1/GP)   ->  6 cols
	bool m_global      = false;	// momento interno global "m" (3x1/GP) ->  6 cols
	// Cinem�tica nos pontos de Gauss
	bool alpha_delta   = false;	// vetor de rota��o no GP (3x1 por GP) ->  6 cols
	bool d_alpha_delta = false;	// taxa de rota��o no GP (3x1 por GP)  ->  6 cols
	bool u_delta       = false;	// deslocamento no GP (3x1 por GP)     ->  6 cols
	bool d_u_delta     = false;	// velocidade no GP (3x1 por GP)       ->  6 cols
	bool A_delta       = false;	// gradiente de deforma��o (3x3/GP)    -> 18 cols
	bool Q_delta       = false;	// tensor de rota��o (3x3 por GP)      -> 18 cols
	bool d_z           = false;	// posi��o deformada diferencial (3x1) ->  6 cols
	bool d_Z           = false;	// tensor anti-sim�trico de d_z (3x3/GP) -> 18 cols
	bool Xi_delta      = false;	// tensor Xi no GP (3x3 por GP)        -> 18 cols
	bool d_A_delta     = false;	// derivada de A_delta (3x3 por GP)    -> 18 cols
	bool d_Xi_delta    = false;	// derivada de Xi_delta (3x3 por GP)   -> 18 cols
	// Fluido interno - Matrix** nos GPs
	bool omega_ip      = false;	// velocidade angular IF (3x1 por GP)  ->  6 cols
	bool temp_if_l     = false;	// carregamento IF resultante (6x1/GP) -> 12 cols
	bool f_pressure_1  = false;	// for�a de press�o 1 (por GP)         -> depende de dim
	bool f_pressure_2  = false;	// for�a de press�o 2 (por GP)         -> depende de dim
	// Din�mica nos GPs
	bool alpha_dot     = false;	// taxa de rota��o din. (3x1 por GP)   ->  6 cols

	// Fun��es de forma nos GPs (Matrix**)
	bool N_gp          = false;	// matriz das fun��es de forma (6x18/GP)  -> 216 cols
	bool deltaN_gp     = false;	// derivadas das fun��es de forma (9x18/GP)-> 324 cols
	// Correnteza mar�tima - Matrix** nos GPs
	bool e3ip          = false;	// dire��o local no GP (3x1 por GP)    ->  6 cols
	bool zi            = false;	// posi��o no GP (3x1 por GP)          ->  6 cols
	bool vel           = false;	// velocidade da correnteza (3x1/GP)   ->  6 cols
	bool velr          = false;	// velocidade relativa (3x1 por GP)    ->  6 cols
	bool element_vel   = false;	// velocidade do elemento (3x1/GP)     ->  6 cols
	bool ut_gp         = false;	// comp. tangencial de velr (3x1/GP)   ->  6 cols
	bool un_gp         = false;	// comp. normal de velr (3x1/GP)       ->  6 cols
	bool d_e3_d_alpha  = false;	// derivada de e3 em rela��o a alpha (3x3/GP) -> 18 cols
	bool Lt            = false;	// operador tangencial (3x3 por GP)    -> 18 cols
	bool Ln            = false;	// operador normal (3x3 por GP)        -> 18 cols
	bool L_u_alpha     = false;	// operador L_u_alpha (3x3 por GP)     -> 18 cols
	bool f_current     = false;	// for�a de correnteza (3x1 por GP)    ->  6 cols
	bool L_gp          = false;	// matriz L (6x6 por GP)               -> 72 cols
	bool t_e           = false;	// vetor tangente unitario (3x1/GP)    ->  6 cols
	bool n_e           = false;	// vetor normal unitario (3x1/GP)      ->  6 cols
	bool vtr           = false;	// velocidade tangencial relativa(3x1/GP)->  6 cols
	bool vnr           = false;	// velocidade normal relativa (3x1/GP) ->  6 cols
	bool Mdt           = false;	// matriz Mdt (3x3 por GP)             -> 18 cols
	bool Mdn           = false;	// matriz Mdn (3x3 por GP)             -> 18 cols
	bool Md2           = false;	// matriz Md2 (6x6 por GP)             -> 72 cols
	// Pipe load - Matrix** nos GPs
	bool kip           = false;	// curvatura atual kip (3x1 por GP)    ->  6 cols
	bool temp_f        = false;	// for�a tempor�ria pipe (3x1/GP)      ->  6 cols
	bool temp_m        = false;	// momento tempor�rio pipe (3x1/GP)    ->  6 cols
	bool temp_l        = false;	// carregamento local pipe (6x1/GP)    -> 12 cols
	bool Kip           = false;	// matriz Kip (3x3 por GP)             -> 18 cols
	bool E3ip          = false;	// matriz E3ip (3x3 por GP)            -> 18 cols
	bool UpsilonN      = false;	// matriz UpsilonN (12x18/GP)          -> 432 cols
	// Din�mica - Matrix** nos GPs
	bool Xi_dot        = false;	// taxa de Xi (3x1 por GP)             ->  6 cols
	bool Mip_gp        = false;	// matriz de massa local (3x3/GP)      -> 18 cols
	bool Jip_gp        = false;	// matriz de in�rcia local (3x3/GP)    -> 18 cols
	bool M_dyn         = false;	// matriz M din�mica (6x6 por GP)      -> 72 cols
	bool Md1           = false;	// matriz Md1 (6x6 por GP)             -> 72 cols

	// --- double* (array nos 2 pontos de Gauss) ---
	bool pressure         = false;	// press�o interna nos GPs             -> 2 cols
	bool rho              = false;	// densidade nos GPs                   -> 2 cols
	bool temperature      = false;	// temperatura nos GPs                 -> 2 cols
	bool flow_velocity    = false;	// velocidade de escoamento nos GPs    -> 2 cols
	bool flow_rate        = false;	// vaz�o nos GPs                       -> 2 cols
	bool internal_pressure= false;	// press�o interna (pipe load) nos GPs -> 2 cols
	bool N1            = false;	// fun��o de forma N1 nos GPs          -> 2 cols
	bool N2            = false;	// fun��o de forma N2 nos GPs          -> 2 cols
	bool N3            = false;	// fun��o de forma N3 nos GPs          -> 2 cols
	bool dN1           = false;	// derivada dN1 nos GPs                -> 2 cols
	bool dN2           = false;	// derivada dN2 nos GPs                -> 2 cols
	bool dN3           = false;	// derivada dN3 nos GPs                -> 2 cols
	bool csi           = false;	// coordenada natural nos GPs          -> 2 cols

	// --- double (escalar) ---
	bool length           = false;	// comprimento indeformado              -> 1 col
	bool jacobian         = false;	// jacobiano                            -> 1 col
	bool alpha1           = false;	// peso de quadratura gaussiana         -> 1 col
	bool p0i              = false;	// press�o interna de refer�ncia (p0i) -> 1 col
	bool p0e              = false;	// press�o externa de refer�ncia (p0e) -> 1 col
	bool rhoi             = false;	// densidade fluido interno             -> 1 col
	bool rhoe             = false;	// densidade fluido externo             -> 1 col
	bool Aint             = false;	// �rea interna da se��o                -> 1 col
	bool rho_adt          = false;	// massa adicional tangencial           -> 1 col
	bool rho_adn          = false;	// massa adicional normal               -> 1 col
	bool load_multiplier  = false;	// multiplicador de carga               -> 1 col
	bool alpha_escalar_delta = false;// escalar de alpha no GP              -> 1 col
	bool g             = false;	// fator g (Rodrigues)                  -> 1 col
	bool signal_t      = false;	// sinal tangencial de Morison          -> 1 col
	bool Cdt           = false;	// coeficiente de arrasto tang.         -> 1 col
	bool Cdn           = false;	// coeficiente de arrasto normal        -> 1 col
	bool Aext          = false;	// �rea externa da se��o                -> 1 col
	bool rho_f         = false;	// densidade do fluido externo          -> 1 col
	bool depth         = false;	// profundidade do elemento             -> 1 col
	bool Un_           = false;	// velocidade normal da correnteza      -> 1 col
	bool un_scalar     = false;	// velocidade normal do elemento        -> 1 col
	bool Ut_           = false;	// velocidade tangencial da correnteza  -> 1 col
	bool ut_scalar     = false;	// velocidade tangencial do elemento    -> 1 col
	bool C1t           = false;	// coeficiente C1 tangencial            -> 1 col
	bool C1n           = false;	// coeficiente C1 normal                -> 1 col
	bool l_factor      = false;	// fator l (steps)                     -> 1 col
	bool mult          = false;	// multiplicador (steps)                -> 1 col
	bool t1            = false;	// contador de step t1                  -> 1 col
	bool t2            = false;	// contador de step t2                  -> 1 col

	// --- Matrix* (matriz �nica) ---
	bool e3r              = false;	// dire��o de refer�ncia local (3x1)   -> 3 cols
	bool e3rg             = false;	// dire��o de refer�ncia global (3x1)  -> 3 cols
	bool omega_rb         = false;	// vel. angular corpo r�gido (3x1)     -> 3 cols
	bool r_inst           = false;	// raio instant�neo                     -> depende de dim
	bool i_loading        = false;	// vetor de esfor�os internos (18x1)   -> 18 cols
	bool e_loading        = false;	// vetor de esfor�os externos (18x1)   -> 18 cols
	bool P_loading        = false;	// vetor de esfor�o desbalanceado(18x1)-> 18 cols
	bool inertial_loading = false;	// esfor�os inerciais (18x1)           -> 18 cols
	bool morison_loading  = false;	// esfor�os de Morison (18x1)          -> 18 cols
	bool transform3       = false;	// matriz de transf. de coord. (3x3)   ->  9 cols
	bool Mr_mat           = false;	// matriz de massa global (dim x dim)  -> depende de dim
	bool Jr_mat           = false;	// matriz de in�rcia global            -> depende de dim
	bool stiffness_mat    = false;	// matriz de rigidez (18x18)           -> 324 cols
	bool D_mat            = false;	// matriz constitutiva (6x6)           -> 36 cols
	bool I3            = false;	// identidade 3x3                       ->  9 cols
	bool B1            = false;	// matriz B1 (6x6)                     -> 36 cols
	bool Qtransp       = false;	// Q transposta (3x3)                  ->  9 cols
	bool B2            = false;	// matriz B2 (6x9)                     -> 54 cols
	bool B2temp        = false;	// matriz B2temp (3x3)                 ->  9 cols
	bool constitutive_stiffness = false; // rigidez constitutiva (18x18) -> 324 cols
	bool geometric_stiffness    = false; // rigidez geom�trica (18x18)   -> 324 cols
	bool loading_stiffness      = false; // rigidez de carregamento(18x18)-> 324 cols
	bool mass          = false;	// matriz de massa (18x18)              -> 324 cols
	bool damping       = false;	// matriz de amortecimento (18x18)      -> 324 cols
	bool mass_modal    = false;	// matriz de massa modal (18x18)        -> 324 cols
	bool damping_modal = false;	// matriz de amortec. modal (18x18)     -> 324 cols
	bool damping_loading= false;	// amortec. de carregamento (18x1)     -> 18 cols
	bool rayleigh_damping= false;	// amortec. de Rayleigh (18x18)        -> 324 cols
	bool transform     = false;	// matriz de transf. completa (18x18)  -> 324 cols
	bool V_alpha_dz_n  = false;	// V(alpha, dz*n) (3x3)                ->  9 cols
	bool V_alpha_m     = false;	// V(alpha, m) (3x3)                   ->  9 cols
	bool d_V_dalpha_apha_m= false;	// dV/dalpha (3x3)                     ->  9 cols
	bool G_d_u_alpha   = false;	// G_d_u_alpha (3x3)                   ->  9 cols
	bool G_d_u_alpha_transp= false;	// G_d_u_alpha transposta (3x3)        ->  9 cols
	bool G_alpha_alpha = false;	// G_alpha_alpha (3x3)                  ->  9 cols
	bool G_alpha_d_alpha= false;	// G_alpha_d_alpha (3x3)               ->  9 cols
	bool G_alpha_d_alpha_transp= false;// G_alpha_d_alpha transposta (3x3) ->  9 cols
	bool B_mat         = false;	// matriz B (6x9)                      -> 54 cols
	bool G_mat         = false;	// matriz G (9x9)                      -> 81 cols
	bool K1ua          = false;	// matriz K1ua (3x3)                   ->  9 cols
	bool K1aa          = false;	// matriz K1aa (3x3)                   ->  9 cols
	bool K2ua          = false;	// matriz K2ua (3x3)                   ->  9 cols
	bool K2au          = false;	// matriz K2au (3x3)                   ->  9 cols
	bool Kext          = false;	// matriz Kext (6x12)                  -> 72 cols
	bool O1            = false;	// matriz O1 (3x3)                     ->  9 cols
	bool Lambda        = false;	// matriz Lambda (3x3)                 ->  9 cols

	// --- vector<double> ---
	bool external_pressure = false;	// press�o externa ao longo do elemento-> N cols (var.)

	// --- LengthPoints (struct: initial + gauss[] + final) ---
	bool length_points     = false;	// posi��es ao longo do elemento: length_points_initial, length_points_gauss_0..N, length_points_final
};

struct UserDefElemMonitorEntry {
	bool          use_element_set = false;	// true = ids s�o ElementSet IDs; false = element IDs
	vector<int>   ids;						// IDs como lidos do input
	UserDefElemMonitorParams params;
	double        time_start = 0.0;			// tempo a partir do qual os dados s�o gravados (opcional)
	vector<int>              resolved_elem_ids;		// IDs de elementos resolvidos - preenchido em StartMonitor()
	vector<UserDefElemMonitorParams> elem_params;	// params por elemento, paralelo a resolved_elem_ids
	vector<double>           elem_time_starts;		// time_start por elemento
	vector<bool>             elem_first_record;		// flag de primeiro registro por elemento
	vector<FILE*>            files;					// um FILE* por elemento resolvido
};
