%% Paper: Thermal Model for Li-Ion Cell %%

%% Geometric parameters

L_p     = 74e-6;   % Length of positive electrode (m)
L_m     = 25e-6;   % Length of membrane separator (m)
L_n     = 75e-6;   % Length of negative electrode (m)
L       = L_p + L_m + L_n;  % Total length (m)

S_p     = 0.0878;    % Geometric surface area of positive electrode (m^2)
S_n     = 0.0878;    % Geometric surface area of negative electrode (m^2)

Rs_p    = 2e-6;     % Radius of particle in positive electrode (m)
Rs_n    = 2e-6;     % Radius of particle in negative electrode (m)

% epse_p  = 0.3;    % Volume fraction of electrolyte in positive electrode (-)
% epse_m0 = 0.45;     % Volume fraction of electrolyte in separator (-)
% epse_n  = 0.4382;    % Volume fraction of electrolyte in negative electrode (-)

epse_p  = 0.338;    % Volume fraction of electrolyte in positive electrode (-)
epse_m0 = 0.37;     % Volume fraction of electrolyte in separator (-)
epse_n  = 0.44;    % Volume fraction of electrolyte in negative electrode (-)

% epsfl_p = 0.15;    % Volume fraction of current conductive fillers in positive electrode (-)
% epsfl_n = 0.0566;     % Volume fraction of current conductive fillers in negative electrode (-)

epsfl_p = 0.142;    % Volume fraction of current conductive fillers in positive electrode (-)
epsfl_n = 0.07;     % Volume fraction of current conductive fillers in negative electrode (-)

epss_p  = (1 - epse_p - epsfl_p);   % Volume fraction of solid active material in positive electrode (-)
epss_n  = (1 - epse_n - epsfl_n);   % Volume fraction of solid active material in negative electrode (-)

as_p    = 3*epss_p/(Rs_p);    % Specific surface area of the porous electrode in positive electrode (m^2/m^3)
as_n    = 3*epss_n/(Rs_n);    % Specific surface area of the porous electrode in negative electrode (m^2/m^3)

%% Material Properties
Ce_avg      = 1000;     % Average concentration of Li ions in electrolyte (mol/m^3)
% Cs_max_p0   = 49943;    % Maximum concentration of solid phase lithium in positive electrode (mol/m^3)
% Cs_max_n0   = 31858;    % Maximum concentration of solid phase lithium in negative electrode (mol/m^3)
Cs_max_p0   = 51555;    % Maximum concentration of solid phase lithium in positive electrode (mol/m^3)
Cs_max_n0   = 30555;    % Maximum concentration of solid phase lithium in negative electrode (mol/m^3)

Ea_dp       = 30e3;      %activation energy for particle diffusion for positive electrode (J/mol);
Ea_dn       = 40e3;      %activation energy for particle diffusion for negative electrode (J/mol);

Ea_rp       = 3.5e4;     %activation energy for reaction for positive electrode (J/mol);
Ea_rn       = 4.5e4;     %activation energy for reaction for negative electrode (J/mol);

Ea_plating  = 3.53e4;
% Ea_plating  = 20e4;
Ea_sei      = 20e4;

Ea_e        = 1.99e4;    %activation energy for electrolyte phase diffusion coefficient
Ea_k        = 1.82e4;    %activation energy for electrolyte phase conductivity

% Ds_p        = (1e-13)*ones(42,1);    % Diffusivity of solid phase in positive electrode (m^2/s)
% Ds_n        = (3.8e-14)*ones(42,1);  % Diffusivity of solid phase in negative electrode (m^2/s)
% De_p        = ones(42,1)*2.5e-10;  % Diffusion coefficient of electrolyte phase in positive electrode (m^2/s)
% De_m        = ones(42,1)*2.5e-10;  % Diffusion coefficient of electrolyte phase in separator (m^2/s)
% De_n        = ones(42,1)*2.5e-10;  % Diffusion coefficient of electrolyte phase in negative electrode (m^2/s)
% kappae_p    = ones(42,1)*2.5;      % Conductivity of electrolyte in positive electrode (S/m)
% kappae_m    = ones(42,1)*2.5;      % Conductivity of electrolyte in separator (S/m)
% kappae_n    = ones(42,1)*2.5;      % Conductivity of electrolyte in negative electrode (S/m)

Ds_p_ref    = 1e-13;    % Diffusion coefficient of solid phase in positive electrode (m^2/s)
Ds_n_ref    = 3.8e-14;  % Diffusion coefficient of solid phase in negative electrode (m^2/s)
% Ds_p_ref    = 1e-11;    % Diffusion coefficient of solid phase in positive electrode (m^2/s)
% Ds_n_ref    = 1.45e-13;  % Diffusion coefficient of solid phase in negative electrode (m^2/s)

De_p_ref    = 2.5e-10;  % Diffusion coefficient of electrolyte phase in positive electrode (m^2/s)
De_m_ref    = 2.5e-10;  % Diffusion coefficient of electrolyte phase in separator (m^2/s)
De_n_ref    = 2.5e-10;  % Diffusion coefficient of electrolyte phase in negative electrode (m^2/s)

kappae_p_ref    = 2.5;      % Conductivity of electrolyte in positive electrode (S/m)
kappae_m_ref    = 2.5;      % Conductivity of electrolyte in separator (S/m)
kappae_n_ref    = 2.5;      % Conductivity of electrolyte in negative electrode (S/m)

sigmas_p    = 10;       % Intrinsic conductivity of active materials in solid phase of positive electrode (1/(Ohm m))
sigmas_n    = 100;      % conductivity of active materials in solid phase of negative electrode (1/(Ohm m))

Brugg_p     = 1.5;      % Fitted parameter
Brugg_m     = 1.5;      % Fitted parameter
Brugg_n     = 1.5;      % Fitted parameter

% Brugg_p     = 1.5;
% Brugg_m     = 2.3;
% Brugg_n     = 4.1;

transf_p      = 0.2*ones(1,42);      % Transference number (-)
transf_m      = 0.2*ones(1,14);
transf_n      = 0.2*ones(1,42);

% transf_p      = 0.435*ones(1,42);      % Transference number (-)
% transf_m      = 0.435*ones(1,14);
% transf_n      = 0.435*ones(1,42);

%% Rate constants

kp_ref          = 0.6e-10;         % Rate constant for positive electrode (m/s)
kn1_ref         = 1.1e-6;        % Rate constant for negative electrode (m/s)
kn2_ref         = 1.1e-12;        % Rate constant for negative electrode (m/s), shif
% kn_plating_ref  = 0;                % Rate constant for Li plating in negative electrode (m/s)
% kn_sei_ref      = 0; 
kn_sei_ref      = 1.1e-16;    % Rate constant for SEI layer formation in negative electrode (m/s)
kn_plating_ref  = 1.1e-17;     % Rate constant for Li plating in negative electrode (m/s)

%% Rate constants

% kp_ref          = 0.6e-10;  % Rate constant for positive electrode (m/s)
% kn1_ref         = 1.1e-6;   % Rate constant for negative electrode (m/s)
% kn2_ref         = 1.1e-12;  % Rate constant for negative electrode (m/s)
% kn_sei_ref      = 20*1.1e-17;    % Rate constant for SEI layer formation in negative electrode (m/s)
% kn_plating_ref  = 200*1.1e-18;    % Rate constant for Li plating in negative electrode (m/s)

%% Constants

F           = 96487;    % Faraday constant (C/mol)
R           = 8.314;    % Universal gas constant (J/K mol)
alphaA      = 0.5;      % Transfer coefficient of an electrochemical rkn_sei_refeaction(-)
rated_curr  = 1.67;     % A 
% I_app           = rated_curr;
%% SEI layer parameters

M_sei       = 0.162;      % Molar mass of SEI (kg/mol)
rho_sei     = 1690;      % Density of SEI (kg/m^3)
k_cond_sei  = 5e-8;  % Conductivity of SEI (S/m)

%% Li plating parameters

M_plating        = 6.94e-3;    % Molar mass of plated Li (kg/mol)
rho_plating      = 0.534e3;      % Density of plated Li (kg/m^3)
k_cond_plating   = 1.1e7;     % Conductivity of plated Li (S/m)

%% Operating conditions

Ta          = 298.15;  	% Temperature (K)
Tref        = 298.15;
EOCV        = 4.2;    	% End of charging voltage (V)
EODV        = 3.0;       % End of discharging voltage (V)
EOCC        = 0.1;       % End of charging current (A)
N_cycles    = 1;        % Number of cycles
% t_cc_char   = 0.01*3600; % Actual time plugged in for cc charging

%% Scale factors
t_sc        = 1; % Scale for x direction
x_sc        = L; % Scale for x direction
Ce_sc       = Ce_avg;   % Scale for electrolyte concentration
Csp_sc      = Cs_max_p0; % Scale for conc of solid phase Li in positive electrode (mol/m^3)
Csn_sc      = Cs_max_n0; % Scale for conc of solid phase Li in negative electrode (mol/m^3)
Csep_sc     = Cs_max_p0; % Scale for conc of solid phase Li in positive electrode - electrolyte boundary (mol/m^3)
Csen_sc     = Cs_max_n0; % Scale for conc of solid phase Li in negative electrode - electrolyte boundary (mol/m^3)
phie_sc     = R*Ta/F;    % Scale for electrolyte potential (V)
phisp_sc    = 1;%I_app * x_sc / (S_p*sigmas_p*epss_p^1.5); % Scale for solid phase potential in positive electrode (V)
phisn_sc    = 1;%I_app * x_sc / (S_n*sigmas_n*epss_n^1.5); % Scale for solid phase potential in positive electrode (V)
% phisp_sc    = I_app * x_sc / (S_p*sigmas_p*epss_p^1.5); % Scale for solid phase potential in positive electrode (V)
% phisn_sc    = I_app * x_sc / (S_n*sigmas_n*epss_n^1.5); % Scale for solid phase potential in positive electrode (V)
%% Solution related params

t_max   = 20*3600;   % maximum time (sec)
N_p     = 21*2;       % Number of nodes in positive electrode 
N_m     = 7*2;       % Number of nodes in separator
N_n     = 21*2;       % Number of nodes in negative electrode 
Nodes   = 1:N_p+N_m+N_n;
Delx_p  = L_p/((N_p - 1)*x_sc);    % Width of node in positive electrode (m)
Delx_m  = L_m/((N_m + 1)*x_sc);    % Width of node in positive electrode (m)
Delx_n  = L_n/((N_n - 1)*x_sc);    % Width of node in positive electrode (m)

%% Energy balance
lambda_m      = 0.16;     %thermal conductivity of the cell 
lambda_p      = 2.1; %thermal conductivity W/m.K
lambda_n      = 1.7; %thermal conductivity W/m.K
% lambda_m      = 1;     %thermal conductivity of the cell 
% lambda_p      = 1; %thermal conductivity W/m.K
% lambda_n      = 1; %thermal conductivity W/m.K
rho_p       = 2328.5; %kg/m^3
rho_n       = 1347.3; %kg/m^3
rho_m       = 1123.0; %kg/m^3
% rho_p       = 2292; %kg/m^3
% rho_n       = 5031.67; %kg/m^3
% rho_p       = 1626; %kg/m^3
% rho_n       = 1626; %kg/m^3
% rho_m       = 1626; %kg/m^3
Cp_p          = 1669.2; %J/kg.K
Cp_m          = 2055.1;
Cp_n          = 1437.4;
% Cp_p          = 750; %J/kg.K
% Cp_m          = 750;
% Cp_n          = 750;
h             = 4.75;%W/m^2.K Natural HT coefficient obt from lit