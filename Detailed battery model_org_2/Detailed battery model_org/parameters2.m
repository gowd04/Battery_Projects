
%% Geometric parameters

L_p     = 7.4e-5;   % Length of positive electrode (m)
L_m     = 2.5e-5;   % Length of membrane separator (m)
L_n     = 7.5e-5;   % Length of negative electrode (m)
L       = L_p + L_m + L_n;  % Total length (m)

S_p     = 0.0878;    % Geometric surface area of positive electrode (m^2)
S_n     = 0.0878;    % Geometric surface area of negative electrode (m^2)

Rs_p    = 2e-6;     % Radius of particle in positive electrode (m)
Rs_n    = 2e-6;     % Radius of particle in negative electrode (m)

epse_p  = 0.338;    % Volume fraction of electrolyte in positive electrode (-)
epse_m0 = 0.37;     % Volume fraction of electrolyte in separator (-)
epse_n  = 0.440;    % Volume fraction of electrolyte in negative electrode (-)

epsfl_p = 0.142;    % Volume fraction of current conductive fillers in positive electrode (-)
epsfl_n = 0.07;     % Volume fraction of current conductive fillers in negative electrode (-)

epss_p  = (1 - epse_p - epsfl_p);   % Volume fraction of solid active material in positive electrode (-)
epss_n  = (1 - epse_n - epsfl_n);   % Volume fraction of solid active material in negative electrode (-)

as_p    = 3*epss_p/(Rs_p);    % Specific surface area of the porous electrode in positive electrode (m^2/m^3)
as_n    = 3*epss_n/(Rs_n);    % Specific surface area of the porous electrode in negative electrode (m^2/m^3)

%% Material Properties
Ce_avg      = 1000;     % Average concentration of Li ions in electrolyte (mol/m^3)
Cs_max_p0   = 51555;    % Maximum concentration of solid phase lithium in positive electrode (mol/m^3)
Cs_max_n0   = 30555;    % Maximum concentration of solid phase lithium in negative electrode (mol/m^3)
Ds_p        = 1e-13;    % Diffusivity of solid phase in positive electrode (m^2/s)
Ds_n        = 3.8e-14;  % Diffusivity of solid phase in negative electrode (m^2/s)
De_p        = 2.5e-10;  % Diffusivity of solid phase in positive electrode (m^2/s)
De_m        = 2.5e-10;  % Diffusivity of solid phase in separator (m^2/s)
De_n        = 2.5e-10;  % Diffusivity of solid phase in negative electrode (m^2/s)
sigmas_p    = 10;       % Intrinsic conductivity of active materials in solid phase of positive electrode (1/(Ohm m))
sigmas_n    = 100;      % Intrinsic conductivity of active materials in solid phase of negative electrode (1/(Ohm m))
kappae_p    = 2.5;      % Conductivity of electrolyte in positive electrode (1/(Ohm m))
kappae_m    = 2.5;      % Conductivity of electrolyte in separator (1/(Ohm m))
kappae_n    = 2.5;      % Conductivity of electrolyte in negative electrode (1/(Ohm m))
Brugg_p     = 1.5;
Brugg_m     = 1.5;
Brugg_n     = 1.5;
transf      = 0.2;      % Transference number (-)

%% Rate constants

kp          = 0.6e-10;  % Rate constant for positive electrode (m/s)
kn1         = 1.1e-6;   % Rate constant for negative electrode (m/s)
kn2         = 1.1e-12;  % Rate constant for negative electrode (m/s)
kn_sei      = 20*1.1e-17;    % Rate constant for SEI layer formation in negative electrode (m/s)
kn_plating  = 200*1.1e-18;    % Rate constant for Li plating in negative electrode (m/s)

%% Constants

F           = 96485;    % Faraday constant (C/mol)
R           = 8.314;    % Universal gas constant (J/K mol)
alphaA      = 0.5;      % Transfer coefficient of an electrochemical reaction(-)
rated_curr  = 1.67;     % A 

%% SEI layer parameters

M_sei       = 0.162;      % Molar mass of SEI (kg/mol)
rho_sei     = 1690;      % Density of SEI (kg/m^3)
k_cond_sei  = 5e-8;  % Conductivity of SEI (S/m)

%% Li plating parameters

M_plating        = 6.94e-3;    % Molar mass of plated Li (kg/mol)
rho_plating      = 0.535e3;      % Density of plated Li (kg/m^3)
k_cond_plating   = 1.1*10^7;     % Conductivity of plated Li (S/m)

%% Operating conditions

T           = 298.15;  	% Temperature (K)
EOCV        = 4.2;    	% End of charging voltage (V)
EODV        = 3.2;        % End of discharging voltage (V)
EOCC        = 0.1;     % End of charging current (A)
N_cycles    = 1;        % Number of cycles
t_cc_char   = 1.1*3600; % Actual time plugged in for cc charging

%% Scale factors

x_sc        = L; % Scale for x direction
t_sc        = 1;  % Scale for time
Ce_sc       = Ce_avg;   % Scale for electrolyte concentration
Csp_sc      = Cs_max_p0; % Scale for conc of solid phase Li in positive electrode (mol/m^3)
Csn_sc      = Cs_max_n0; % Scale for conc of solid phase Li in negative electrode (mol/m^3)
Csep_sc     = Cs_max_p0; % Scale for conc of solid phase Li in positive electrode - electrolyte boundary (mol/m^3)
Csen_sc     = Cs_max_n0; % Scale for conc of solid phase Li in negative electrode - electrolyte boundary (mol/m^3)
phie_sc     = R*T/F;    % Scale for electrolyte potential (V)
phisp_sc    = 1;%I_app * x_sc / (S_p*sigmas_p*epss_p^1.5); % Scale for solid phase potential in positive electrode (V)
phisn_sc    = 1;%I_app * x_sc / (S_n*sigmas_n*epss_n^1.5); % Scale for solid phase potential in positive electrode (V)

%% Solution related params

t_max   = 5*3600;   % maximum time (sec)
N_p     = 21*2;       % Number of nodes in positive electrode 
N_m     = 7*2;       % Number of nodes in separator
N_n     = 21*2;       % Number of nodes in negative electrode 
Delx_p  = L_p/((N_p - 1)*x_sc);    % Width of node in positive electrode (m)
Delx_m  = L_m/((N_m + 1)*x_sc);    % Width of node in positive electrode (m)
Delx_n  = L_n/((N_n - 1)*x_sc);    % Width of node in positive electrode (m)



