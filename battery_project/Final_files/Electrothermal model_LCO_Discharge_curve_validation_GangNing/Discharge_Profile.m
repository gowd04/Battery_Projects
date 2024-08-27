%% Paper: Thermal Model for Li-Ion Cell %%

clc
clear
% close all;

tic
parameters3
initial_condns_discharge

load('possible_sites.mat');
data = csvread('disch_curve_gang_ning.csv');

%% Initializations

time = 0;
C = [0.5];
% Temp_cell          = [];
% Vcell              = [];
% Current_cell       = [];
% time_tot           = [];
for i = 1:length(C)
% Temp_cell           = [];
% Vcell               = [];
% Current             = [];
% Ce                  = [];
% Pe                  = [];
% Csp                 = [];
% Csn                 = [];
% Csep                = [];
% Csen                = [];
% Psp                 = [];
% Psn                 = [];
% del_sei_all         = [0 delta_sei*ones(1,N_n)];
% R_sei_all           = [0 R_sei*ones(1,N_n)];
% C_loss_sei_all      = [0 C_loss_sei];
% C_loss_plating_all  = [0 C_loss_plating];
% R_plating_all       = [0 R_plating*ones(1,N_n)];

%% Discretization in x
xp       = 0:Delx_p:L_p/x_sc;
xm       = (L_p/x_sc):Delx_m:(L_p+L_m)/x_sc;
xn       = ((L_p + L_m)/x_sc):Delx_n:(L_p + L_m+L_n)/x_sc;

x_scaled = [xp'; xm(2:end-1)';xn'];

%% Initial condition

I_app       = 0.5*rated_curr;
% I_app       = C(i)*rated_curr;
C1_p        = Cs0_p*ones(1,N_p)/Csp_sc;
C1_n        = Cs0_n*ones(1,N_n)/Csn_sc;
C2_p        = Ce0*ones(1,N_p)/Ce_sc; %electrolyte phase conc of Li
C2_m        = Ce0*ones(1,N_m)/Ce_sc; %electrolyte phase conc of Li
C2_n        = Ce0*ones(1,N_n)/Ce_sc; %electrolyte phase conc of Li

T_p         = T0*(ones(1,N_p))/Ta;
T_s         = T0*(ones(1,N_m))/Ta;
T_n         = T0*(ones(1,N_n))/Ta;

Jptot       = I_app/(as_p*S_p*L_p);
Jntot       = -I_app/(as_n*S_n*L_n);

C3_p        = C1_p*Csp_sc/Csep_sc -(Jptot*Rs_p/(5*F*Ds_p_ref*Csep_sc)); %solid phase conc at interphase
C3_n        = C1_n*Csn_sc/Csen_sc -(Jntot*Rs_n/(5*F*Ds_n_ref*Csen_sc));

Cs_max_p    = Cs_max_p0 - C_loss_sei - C_loss_plating;
Cs_max_n    = Cs_max_n0 - C_loss_sei - C_loss_plating;

theta_p0    = C3_p*Csep_sc/Cs_max_p;
theta_n0    = C3_n*Csen_sc/Cs_max_n;

Uref0_p     = 1.654107793108376e+06 *theta_p0.^10 -1.249511527828579e+07 *theta_p0.^9 + 4.215812681219930e+07*theta_p0.^8 -8.365902527257702e+07*theta_p0.^7 + ...
    1.081256432503016e+08*theta_p0.^6 -9.510080800780240e+07*theta_p0.^5 + 5.764438732509381e+07*theta_p0.^4 -2.377606148390258e+07*theta_p0.^3 + ...
    6.386329451323811e+06*theta_p0.^2 -1.008737242075382e+06*theta_p0 +  7.115616223683314e+04;

Uref0_n     = (1.19970 + (118.1911 * (theta_n0.^0.5)) - (706.0711 * (theta_n0)) + (2217.6479 * (theta_n0.^1.5)) - (1675.1321 * (theta_n0.^2)))./...
    (1.0 + (131.7572 * (theta_n0.^0.5)) - (32.1402 * (theta_n0)) - (746.8463 * (theta_n0.^1.5)) + (15502.9505 * (theta_n0.^2)) - (14213.0747 * (theta_n0.^2.5)));

% Uref0_p = (-4.656 + (88.669*(theta_p0.^2)) - (401.119*(theta_p0.^4)) + (342.909*(theta_p0.^6)) - (462.471*(theta_p0.^8)) + (433.434*(theta_p0.^10)))./...
%     (-1 + (18.933*(theta_p0.^2)) - (79.532*(theta_p0.^4)) + 37.311*(theta_p0.^6) - (73.083*(theta_p0.^8)) + (95.96*(theta_p0.^10)));
% 
% Uref0_n = 0.7222 + 0.1387*theta_n0 + 0.029*(theta_n0.^0.5) - (0.0172./theta_n0) + (0.0019./(theta_n0.^1.5)) + 0.2808*exp(0.9-(15*theta_n0))...
%     -0.7984*exp((0.4465*theta_n0) - 0.4108);

P1p0        = Uref0_p/phisp_sc;%
P1n0        = Uref0_n/phisn_sc;%
P20         = zeros(1,N_p+N_m+N_n);

Jpi_init    = ones(1,N_p);
Jni_init    = ones(1,N_n);

epse_m      = (epse_m0*S_p*L_m - sum(Li_den_matrix(:,5)))/(S_p*L_m);
epse_m_all  = [0 epse_m];

y           = [C1_p(end,:), C1_n(end,:), C2_p(end,:),C2_m(end,:),C2_n(end,:),P1p0, P1n0, P20,C3_p(end,:),C3_n(end,:), T_p(end,:), T_s(end,:), T_n(end,:) Jpi_init, Jni_init];

%% CC discharging (Index_chdisch = 1)
Index_chdisch = 1;

I_app   = -C(i)*rated_curr;
% t_max   = t_max*10;
y0_all  = y(end,:);
options = odeset('Mass',@massmatrix_gen_DAE,'RelTol',1e-2,'MaxStep',100,'Events',@breakode_battery_new);
[t,y]   = ode15s(@battery_pde_model_DAE_new2,[0 t_max/t_sc],y0_all,options,Index_chdisch,I_app,[],R_sei, R_plating, Cs_max_p,Cs_max_n,epse_m);

time = t*t_sc;

C1_p = y(:,1:N_p);
C1_n = y(:,N_p+1:N_p+N_n);
C2_p = y(:,N_p+N_n+1:2*N_p+N_n);
C2_m = y(:,2*N_p+N_n+1:2*N_p+N_m+N_n);
C2_n = y(:,2*N_p+N_m+N_n+1:2*N_p+N_m+2*N_n);
P1_p = y(:,2*N_p+N_m+2*N_n+1:3*N_p+N_m+2*N_n);
P1_n = y(:,3*N_p+N_m+2*N_n+1:3*N_p+N_m+3*N_n);
P2_p = y(:,3*N_p+N_m+3*N_n+1:4*N_p+N_m+3*N_n);
P2_m = y(:,4*N_p+N_m+3*N_n+1:4*N_p+2*N_m+3*N_n);
P2_n = y(:,4*N_p+2*N_m+3*N_n+1:4*N_p+2*N_m+4*N_n);
C3_p = y(:,4*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+4*N_n);
C3_n = y(:,5*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+5*N_n);
T_p  = y(:,5*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+5*N_n);
T_s  = y(:,6*N_p+2*N_m+5*N_n+1:6*N_p+3*N_m+5*N_n);
T_n  = y(:,6*N_p+3*N_m+5*N_n+1:6*N_p+3*N_m+6*N_n);

% kp = kp_ref*exp((-Ea_rp/R)*((1./(T_p*Ta))-(1/Ta))); % rate constant for positive electrode
% kn1 = kn1_ref*exp((-Ea_rn/R)*((1./(T_n*Ta))-(1/Ta)));
% kn2 = kn2_ref*exp((-Ea_rp/R)*((1./(T_p*Ta))-(1/Ta)));

kp  = kp_ref*exp((-Ea_rp/R)*((1./(Ta))-(1/Tref))); % rate constant for positive electrode
kn1 = kn1_ref*exp((-Ea_rn/R)*((1./(Ta))-(1/Tref)));
kn2 = kn2_ref*exp((-Ea_rp/R)*((1./(Ta))-(1/Tref)));
kn_sei  = kn_sei_ref*exp((-Ea_sei/R)*((1./(Ta))-(1/Tref)));
kn_plating = kn_plating_ref*exp((-Ea_rn/R)*((1./(Ta))-(1/Tref)));
% kp = mean(kp_arr);
% kn1 = mean(kn1_arr);
% kn2 = mean(kn2_arr);
% kn_sei = mean(kn_sei_arr);
% kn_plating = mean(kn_plating_arr);

Jpi0      = (kp * F *((C3_p*Csep_sc).^(1-alphaA)).*((C2_p*Ce_sc).^alphaA));
Jni0_des  = (((kn1)^(1-alphaA))* ((kn2)^(1-alphaA)) * F*((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));

Jpi       = y(:,6*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+6*N_n)*abs(I_app/S_p);
Jni_des   = y(:,7*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+7*N_n)*abs(I_app/S_n);

Vcell   = P1_p(:,1)*phisp_sc - P1_n(:,end)*phisn_sc;
Current = I_app*ones(length(t),1);
Ce      = [C2_p, C2_m,C2_n]*Ce_sc;
Pe      = [P2_p, P2_m,P2_n]*phie_sc;
Csp     = C1_p*Csp_sc;
Csn     = C1_n*Csn_sc;
Csep    = C3_p*Csep_sc;
Csen    = C3_n*Csen_sc;
Psp     = P1_p*phisp_sc;
Psn     = P1_n*phisn_sc;
Temp   = [T_p T_s T_n]*T0;
% Temp_cell = [Temp_cell Temp];
% Current_cell = [Current_cell Current];
% time_tot = [time_tot time];
Simulation_time = toc
% clearvars -except Simulation_time Ce Pe Csp Csn Csep Csen Psp Psn x_scaled Temp Nodes N_cycles time Vcell Current C_loss_sei_all del_sei_all R_sei_all Li_den_matrix R_plating_all C_loss_plating_all epse_m_all
% clearvars -except Temp data data2 Nodes N_cycles time Vcell Current C_loss_sei_all del_sei_all R_sei_all Li_den_matrix R_plating_all C_loss_plating_all epse_m_all
Simulation_time
%% Generating plots

% figure
% plot(Nodes, Temp_cell(end,:))
% xlabel('Node')
% ylabel('Temperature')

% figure
% plot(time/3600, Temp_cell(:,1))
% xlabel('Node')
% ylabel('Temperature')

%Temperature with discharge capacity
% subplot(2,1,1);
% hold on
% plot(data(:,10),data(:,11),'or',((time.*abs(Current))/(3600)), Temp(:,1));
% 
% xlabel('Discharge Capacity (Ah)')
% ylabel('Temperature at x=0')
% % legend('1C','simulation', 'Location', 'NorthEast')
% ylim([298 302])

% Cell voltage with discharge capacity
% subplot(2,1,2);
% hold on
figure
plot(data(:,1),data(:,2),'or',((time.*abs(Current))/(3600.0)) , Vcell);
xlabel('Discharge capacity (Ah)')
ylabel('Cell voltage (V)')
% legend('1C','simulation', 'Location', 'NorthEast')
ylim([3.2 4.3])

% hold on
% plot(time.*abs(Current)/(3600), Temp_cell(:,end))
% xlabel('Discharge capacity (Ah)')
% ylabel('Temperature at x=L')
% hold off
%     figure
%     plot(data2(:,1),data2(:,2),'om',time, Vcell)
%     xlabel('Time (s)')
%     ylabel('Cell voltage (V)')

% % Electrolyte concentration with distance at various times
% figure
% plot(x_scaled,Ce)
% xlabel('Dimensionless distance (x)')
% ylabel('Concentration (mol/m^3)')
% title('Electrolyte concentration with distance at various times')
% 
% % Electrolyte potential with distance at various times
% figure
% plot(x_scaled,Pe)
% xlabel('Dimensionless distance (x)')
% ylabel('Potential (V)')
% title('Electrolyte potential with distance at various times')
% 
% % Solid phase concentration in positive electrode
% figure
% plot(xp,Csp)
% xlabel('Dimensionless distance (x)')
% ylabel('Concentration (mol/m^3)')
% title('C_{s,p} variation with distance at various times')
% 
% % Solid phase concentration in negative electrode
% figure
% plot(xn - xn(1),Csn)
% xlabel('Dimensionless distance (x)')
% ylabel('Concentration (mol/m^3)')
% title('C_{s,n} variation with distance at various times')
% 
% % Electrolyte concentration variation with time
% figure
% plot(time/3600,Ce)
% xlabel('Time, sec')
% ylabel('Concentration (mol/m^3)')
% title('Electrolyte concentration with time at various positions')
% 
% % Electrolyte potential variation with time
% figure
% plot(time/3600,Pe)
% xlabel('Time, sec')
% ylabel('Potential (V)')
% title('Electrolyte potential with time at various positions')

end