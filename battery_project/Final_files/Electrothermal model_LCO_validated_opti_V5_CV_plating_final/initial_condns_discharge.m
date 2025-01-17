%% Initial conditions

Ce0            = Ce_avg;              % Initial concentration of Li ions in electrolyte (mol/m^3)
T0 = Ta; % initial condition for temperature(K)
rated_curr  = 1.67;     % A
%U0 = (P1p0 + P1n0)/2;
Cs0_p          = 0.48*Cs_max_p0;    % Initial concentration of solid phase Li in positive electrode (mol/m^3)
Cs0_n          = 0.83*Cs_max_n0;
delta_sei      = 0;        % Inital SEI film thickness(Ohm/m2)
R_sei          = 0*L_m;    % Inital resistance due to SEI film (Ohm*m2)
C_loss_sei     = 0;

R_plating      = 0;
C_loss_plating = 0;
