clc;
clear;
% close all
tic
% cd E:\Research\Codes\ElectroThermal battery model
parameters3
initial_condns_discharge
load('possible_sites.mat')

%% Initializations

time                = 0;
Temp_cell           = [];
Temp_avg            = [];
Vcell               = [];
Current             = [];
Ce                  = [];
Pe                  = [];
Csp                 = [];
Csn                 = [];
Csep                = [];
Csen                = [];
Psp                 = [];
Psn                 = [];
del_sei_all         = [0 delta_sei*ones(1,N_n)];
R_sei_all           = [0 R_sei*ones(1,N_n)];
C_loss_sei_all      = [0 C_loss_sei];      
C_loss_plating_all  = [0 C_loss_plating];
R_plating_all       = [0 R_plating*ones(1,N_n)];

%% Discretization in x
xp       = 0:Delx_p:L_p/x_sc;
xm       = (L_p/x_sc):Delx_m:(L_p+L_m)/x_sc;
xn       = ((L_p + L_m)/x_sc):Delx_n:(L_p + L_m+L_n)/x_sc;
x_scaled = [xp'; xm(2:end-1)';xn'];

%% Initial condition 

I_app       = 1.0*rated_curr;
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

P1p0        = Uref0_p/phisp_sc;
P1n0        = Uref0_n/phisn_sc;
P20         = zeros(1,N_p+N_m+N_n);
% U0 = (P1p0 + P1n0)/2; %open circuit voltage
Jpi_init    = ones(1,N_p);
Jni_init    = ones(1,N_n);

epse_m      = (epse_m0*S_p*L_m - sum(Li_den_matrix(:,5)))/(S_p*L_m);
epse_m_all  = [0 epse_m];

y           = [C1_p(end,:), C1_n(end,:), C2_p(end,:),C2_m(end,:),C2_n(end,:),P1p0, P1n0, P20,C3_p(end,:),C3_n(end,:), T_p(end,:), T_s(end,:), T_n(end,:), Jpi_init, Jni_init];

%% Charging discharging cycles

for ndx = 1:N_cycles
    disp(['Cycle No. = ',num2str(ndx)]) 
    
    %% CC charging (Index_chdisch = 2)
    
%     Index_chdisch   = 2;
%     I_app           = 1.5*rated_curr;
%     y0_all          = [y(end,:), y(end,end-N_n+1:end)] %changed the index to not include the columns belonging to T variable
%     options         = odeset('Mass',@massmatrix_gen_DAE,'MaxStep',100,'RelTol',1e-2,'Events',@breakode_battery);
%     [t,y]           = ode23t(@battery_pde_model_DAE_new,[0 t_max/t_sc],y0_all,options,Index_chdisch, I_app,Vcell, R_sei, R_plating, Cs_max_p,Cs_max_n,epse_m);
    
%     time            = [time; time(end) + (t*t_sc)];
%   
%     C1_p            = y(:,1:N_p);
%     C1_n            = y(:,N_p+1:N_p+N_n);
%     C2_p            = y(:,N_p+N_n+1:2*N_p+N_n);
%     C2_m            = y(:,2*N_p+N_n+1:2*N_p+N_m+N_n);
%     C2_n            = y(:,2*N_p+N_m+N_n+1:2*N_p+N_m+2*N_n);
%     P1_p            = y(:,2*N_p+N_m+2*N_n+1:3*N_p+N_m+2*N_n);
%     P1_n            = y(:,3*N_p+N_m+2*N_n+1:3*N_p+N_m+3*N_n);
%     P2_p            = y(:,3*N_p+N_m+3*N_n+1:4*N_p+N_m+3*N_n);
%     P2_m            = y(:,4*N_p+N_m+3*N_n+1:4*N_p+2*N_m+3*N_n);
%     P2_n            = y(:,4*N_p+2*N_m+3*N_n+1:4*N_p+2*N_m+4*N_n);
%     C3_p            = y(:,4*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+4*N_n);
%     C3_n            = y(:,5*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+5*N_n);
%     T_p             = y(:,5*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+5*N_n);
%     T_s             = y(:,6*N_p+2*N_m+5*N_n+1:6*N_p+3*N_m+5*N_n);
%     T_n             = y(:,6*N_p+3*N_m+5*N_n+1:6*N_p+3*N_m+6*N_n);
%     
%     Jni_des         = y(:,7*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+7*N_n)*abs(I_app/S_n);
%     Jni_sei         = y(:,7*N_p+3*N_m+7*N_n+1:7*N_p+3*N_m+8*N_n)*abs(I_app/S_n);
%     
%     delta_sei       = sum(abs(Jni_sei-Jni_des)*M_sei.*repmat([0;diff(t*t_sc)],1,N_n)/(rho_sei*F)) + delta_sei;
%     R_sei           = delta_sei/k_cond_sei;
%     
%     C_loss_sei      = sum(sum(abs(Jni_sei-Jni_des).*repmat([0;diff(t*t_sc)],1,N_n)*as_n/(F))) + C_loss_sei;
%     
%     Cs_max_p        = Cs_max_p0 - C_loss_sei;
%     Cs_max_n        = Cs_max_n0 - C_loss_sei;
%     
%     Vcell           = [Vcell; P1_p(:,1)*phisp_sc - P1_n(:,end)*phisn_sc];
%     Current         = [Current; I_app*ones(length(t),1)];
%     Ce              = [Ce; [C2_p, C2_m,C2_n]*Ce_sc];
%     Pe              = [Pe; [P2_p, P2_m,P2_n]*phie_sc];
%     Csp             = [Csp; C1_p*Csp_sc];
%     Csn             = [Csn; C1_n*Csn_sc];
%     Csep            = [Csep; C3_p*Csep_sc];
%     Csen            = [Csen; C3_n*Csen_sc];
%     Psp             = [Psp; P1_p*phisp_sc];
%     Psn             = [Psn; P1_n*phisn_sc];
%     C_loss_sei_all  = [C_loss_sei_all; time(end) C_loss_sei];
%     del_sei_all     = [del_sei_all; time(end) delta_sei];
%     R_sei_all       = [R_sei_all; time(end) R_sei];
%     Temp_cell       = [T_p T_s T_n]*Ta;
%     Tcell           = 
    
%% CC charging with Li plating (Index_chdisch = 4)
    
        Index_chdisch   = 4;
        y0_all          = [y(end,:), y(end,end-(N_n+N_p-1):end)]; %changed the index to not include the columns belonging to T variable
        options         = odeset('Mass',@massmatrix_gen_DAE,'MaxStep',100,'RelTol',1e-2,'Events',@breakode_battery_new);
        [t,y]           = ode23t(@battery_pde_model_DAE_new2,[0 t_max/t_sc], y0_all,options,Index_chdisch, I_app,[], R_sei, R_plating, Cs_max_p,Cs_max_n,epse_m);
        time            = [time; time(end) + (t*t_sc)];
        C1_p            = y(:,1:N_p);
        C1_n            = y(:,N_p+1:N_p+N_n);
        C2_p            = y(:,N_p+N_n+1:2*N_p+N_n);
        C2_m            = y(:,2*N_p+N_n+1:2*N_p+N_m+N_n);
        C2_n            = y(:,2*N_p+N_m+N_n+1:2*N_p+N_m+2*N_n);
        P1_p            = y(:,2*N_p+N_m+2*N_n+1:3*N_p+N_m+2*N_n);
        P1_n            = y(:,3*N_p+N_m+2*N_n+1:3*N_p+N_m+3*N_n);
        P2_p            = y(:,3*N_p+N_m+3*N_n+1:4*N_p+N_m+3*N_n);
        P2_m            = y(:,4*N_p+N_m+3*N_n+1:4*N_p+2*N_m+3*N_n);
        P2_n            = y(:,4*N_p+2*N_m+3*N_n+1:4*N_p+2*N_m+4*N_n);
        C3_p            = y(:,4*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+4*N_n);
        C3_n            = y(:,5*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+5*N_n);
        T_p             = y(:,5*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+5*N_n);
        T_s             = y(:,6*N_p+2*N_m+5*N_n+1:6*N_p+3*N_m+5*N_n);
        T_n             = y(:,6*N_p+3*N_m+5*N_n+1:6*N_p+3*N_m+6*N_n);
        
        Jni_des         = y(:,7*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+7*N_n)*abs(I_app/S_n);
        Jni_sei         = y(:,7*N_p+3*N_m+7*N_n+1:7*N_p+3*N_m+8*N_n)*abs(I_app/S_n);
        Jni_sei_plating = y(:,7*N_p+3*N_m+8*N_n+1:7*N_p+3*N_m+9*N_n)*abs(I_app/S_n);
               
        delta_sei       = sum(abs(Jni_sei-Jni_des)*M_sei.*repmat([0;diff(t*t_sc)],1,size(Jni_sei,2))/(rho_sei*F)) + delta_sei;
        R_sei           = delta_sei/k_cond_sei;
        
        C_loss_sei      = sum(sum(abs(Jni_sei-Jni_des).*repmat([0;diff(t*t_sc)],1,N_n)*as_n/F)) + C_loss_sei;
        C_loss_plating  = sum(sum(abs(Jni_sei_plating-Jni_sei).*repmat([0;diff(t*t_sc)],1,N_n)*as_n/F)) + C_loss_plating;
        
        Li_plated_vol   = sum(sum(abs(Jni_sei_plating-Jni_sei).*repmat([0;diff(t*t_sc)],1,N_n)*as_n*S_n*L_n*M_plating/(F*rho_plating)));
        
        Li_den_matrix   = Lithium_plating_dendrites_new(Li_den_matrix, Li_plated_vol);
        
        dendrite_vol    = sum(Li_den_matrix(:,5));

        % Volume fraction of electrolyte in separator decreases as the
        % dendritic volume increases. This is represented by the eq. below

        epse_m          = (epse_m0*S_p*L_m - dendrite_vol)/(S_p*L_m);
        
        R_plating(1)    = sum(Li_den_matrix(:,4)./(k_cond_plating));
        R_plating(2:N_n)= 0;
        
        Cs_max_p        = Cs_max_p0 - C_loss_sei - C_loss_plating;
        Cs_max_n        = Cs_max_n0 - C_loss_sei - C_loss_plating;
        
        Vcell           = [Vcell; P1_p(:,1)*phisp_sc - P1_n(:,end)*phisn_sc];
        Current         = [Current; I_app*ones(length(t),1)];
        Ce              = [Ce; [C2_p, C2_m,C2_n]*Ce_sc];
        Pe              = [Pe; [P2_p, P2_m,P2_n]*phie_sc];
        Csp             = [Csp; C1_p*Csp_sc];
        Csn             = [Csn; C1_n*Csn_sc];
        Csep            = [Csep; C3_p*Csep_sc];
        Csen            = [Csen; C3_n*Csen_sc];
        Psp             = [Psp; P1_p*phisp_sc];
        Psn             = [Psn; P1_n*phisn_sc];
        del_sei_all     = [del_sei_all; time(end) delta_sei];
        R_sei_all       = [R_sei_all; time(end) R_sei];
        Temp_cell       = [Temp_cell;[T_p T_s T_n]*T0];
        
        C_loss_sei_all      = [C_loss_sei_all; time(end) C_loss_sei];
        C_loss_plating_all  = [C_loss_plating_all; time(end) C_loss_plating];
        R_plating_all       = [R_plating_all; time(end) R_plating];
        epse_m_all          = [epse_m_all; time(end) epse_m];
        
        y                   = y(:,1:end-N_n);
%     end

    %% CV charging (Index_chdisch = 3)
    
    Index_chdisch = 3;
    y0_all  = [y(end,:), I_app];
    
    options = odeset('Mass',@massmatrix_gen_DAE,'RelTol',1e-2,'MaxStep',100,'Events',@breakode_battery_new);
    [t,y]   = ode23t(@battery_pde_model_DAE_new2,[0 0.25*t_max/t_sc],y0_all,options,Index_chdisch,[],Vcell(end),R_sei, R_plating, Cs_max_p,Cs_max_n,epse_m);
    
    time = [time; time(end) + (t*t_sc)];
    
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
    %     Jni_des   = y(:,6*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+6*N_n)*abs(I_app/S_n);
    %     Jni_sei   = y(:,6*N_p+2*N_m+6*N_n+1:6*N_p+2*N_m+7*N_n)*abs(I_app/S_n);
    I_app = y(:,7*N_p+3*N_m+8*N_n+1);
    
    Jni_des   = y(:,7*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+7*N_n).*abs(repmat(I_app,1,N_n)/S_n);
    Jni_sei   = y(:,7*N_p+3*N_m+7*N_n+1:7*N_p+3*N_m+8*N_n).*abs(repmat(I_app,1,N_n)/S_n);
          
    delta_sei  	= sum(abs(Jni_sei-Jni_des)*M_sei.*repmat([0;diff(t*t_sc)],1,N_n)/(rho_sei*F)) + delta_sei;
    R_sei       = delta_sei/k_cond_sei;
    
    C_loss_sei  = sum(sum(abs(Jni_sei-Jni_des).*repmat([0;diff(t*t_sc)],1,N_n)*as_n/(F))) + C_loss_sei;
    
    Cs_max_p    = Cs_max_p0 - C_loss_sei - C_loss_plating;
    Cs_max_n    = Cs_max_n0 - C_loss_sei - C_loss_plating;
    
    Vcell       = [Vcell; P1_p(:,1)*phisp_sc - P1_n(:,end)*phisn_sc];
    Current     = [Current; y(:,6*N_p+2*N_m+7*N_n+1)];
    Ce          = [Ce; [C2_p, C2_m,C2_n]*Ce_sc];
    Pe          = [Pe; [P2_p, P2_m,P2_n]*phie_sc];
    Csp         = [Csp; C1_p*Csp_sc];
    Csn         = [Csn; C1_n*Csn_sc];
    Csep        = [Csep; C3_p*Csep_sc];
    Csen        = [Csen; C3_n*Csen_sc];
    Psp         = [Psp; P1_p*phisp_sc];
    Psn         = [Psn; P1_n*phisn_sc];
    C_loss_sei_all  = [C_loss_sei_all; time(end) C_loss_sei];
    del_sei_all = [del_sei_all; time(end) delta_sei];
    R_sei_all   = [R_sei_all; time(end) R_sei];    
    Temp_cell   = [Temp_cell;[T_p T_s T_n]*T0];
    Temp_avg    = mean(Temp_cell');
    Temp        = Temp_avg';
    
    %% CC discharging (Index_chdisch = 1)    
    
%     Index_chdisch = 1;
%     
%     I_app   = -0.5*rated_curr;
%     y0_all  = y(end,1:end-1-N_n);
%     
%     options = odeset('Mass',@massmatrix_gen_DAE,'RelTol',1e-2,'MaxStep',100,'Events',@breakode_battery_new);
%     [t,y]   = ode15s(@battery_pde_model_DAE_new2,[0 t_max/t_sc],y0_all,options,Index_chdisch,I_app,[],R_sei, R_plating, Cs_max_p,Cs_max_n,epse_m);
%     
%     time = [time; time(end) + (t*t_sc)];
%     
%     C1_p = y(:,1:N_p);
%     C1_n = y(:,N_p+1:N_p+N_n);
%     C2_p = y(:,N_p+N_n+1:2*N_p+N_n);
%     C2_m = y(:,2*N_p+N_n+1:2*N_p+N_m+N_n);
%     C2_n = y(:,2*N_p+N_m+N_n+1:2*N_p+N_m+2*N_n);
%     P1_p = y(:,2*N_p+N_m+2*N_n+1:3*N_p+N_m+2*N_n);
%     P1_n = y(:,3*N_p+N_m+2*N_n+1:3*N_p+N_m+3*N_n);
%     P2_p = y(:,3*N_p+N_m+3*N_n+1:4*N_p+N_m+3*N_n);
%     P2_m = y(:,4*N_p+N_m+3*N_n+1:4*N_p+2*N_m+3*N_n);
%     P2_n = y(:,4*N_p+2*N_m+3*N_n+1:4*N_p+2*N_m+4*N_n);
%     C3_p = y(:,4*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+4*N_n);
%     C3_n = y(:,5*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+5*N_n);
%     T_p  = y(:,5*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+5*N_n);
%     T_s  = y(:,6*N_p+2*N_m+5*N_n+1:6*N_p+3*N_m+5*N_n);
%     T_n  = y(:,6*N_p+3*N_m+5*N_n+1:6*N_p+3*N_m+6*N_n);
%     
%     Jpi0      = (kp_ref * F*((C3_p*Csep_sc).^(1-alphaA)).*((C2_p*Ce_sc).^alphaA));
%     Jni0_des  = (((kn1_ref)^(1-alphaA))* ((kn2_ref)^(1-alphaA)) * F*((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));
%     Jpi       = y(:,6*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+6*N_n)*abs(I_app/S_p);
%     Jni_des   = y(:,7*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+7*N_n)*abs(I_app/S_n);
%         
%     Vcell   = [Vcell; P1_p(:,1)*phisp_sc - P1_n(:,end)*phisn_sc];
%     Current = [Current; I_app*ones(length(t),1)];
%     Ce      = [Ce; [C2_p, C2_m,C2_n]*Ce_sc];
%     Pe      = [Pe; [P2_p, P2_m,P2_n]*phie_sc];
%     Csp     = [Csp; C1_p*Csp_sc];
%     Csn     = [Csn; C1_n*Csn_sc];
%     Csep    = [Csep; C3_p*Csep_sc];
%     Csen    = [Csen; C3_n*Csen_sc];
%     Psp     = [Psp; P1_p*phisp_sc];
%     Psn     = [Psn; P1_n*phisn_sc];
%     Temp_cell   = [Temp_cell;[T_p T_s T_n]*T0];
%     Temp_avg    = mean(Temp_cell');
%     Temp        = Temp_avg';
    
end
time = time(2:end);
Simulation_time = toc;
% clearvars -except Simulation_time Ce Pe Csp Csn Csep Csen Psp Psn x_scaled Temp_cell Nodes N_cycles time Vcell Current C_loss_sei_all del_sei_all R_sei_all Li_den_matrix R_plating_all C_loss_plating_all epse_m_all
% save(['Results_',num2str(N_cycles),'cycles_without_SEI_plating.mat'])
Simulation_time
%%
% figure
% plot(Psp(:,end))
% 
% figure
% plot(Psn(:,1))
% 
% figure
% plot(Vcell)

%% Generating plots

% figure
% plot(Nodes, Temp_cell(end,:))
% xlabel('Node')
% ylabel('Temperature').x

figure
plot(time/3600, Temp_avg)
xlabel('Time')
ylabel('Temperature')

% Cell voltage with discharge capacity
% figure
% plot(time.*abs(Current)/3600, Vcell)
% xlabel('Discharge capacity (Ah)')
% ylabel('Cell voltage (V)')

% Cell voltage with time
figure
plot(time/3600, Vcell)
xlabel('Time (hr)')
ylabel('Cell voltage (V)')

% % Electrolyte concentration with distance at various times
% figure
% plot(x_scaled,[C2_p, C2_m,C2_n]'*Ce_sc)
% xlabel('Dimensionless distance (x)')
% ylabel('Concentration (mol/m^3)')
% title('Electrolyte concentration with distance at various times')
% 
% % Electrolyte potential with distance at various times
% figure
% plot(x_scaled,[P2_p,P2_m,P2_n]'*phie_sc)
% xlabel('Dimensionless distance (x)')
% ylabel('Potential (V)')
% title('Electrolyte potential with distance at various times')
% 
% % Solid phase concentration in positive electrode
% figure
% plot(xp,C1_p*Csp_sc)
% xlabel('Dimensionless distance (x)')
% ylabel('Concentration (mol/m^3)')
% title('C_{s,p} variation with distance at various times')
% 
% % Solid phase concentration in negative electrode
% figure
% plot(xn - xn(1),C1_n*Csn_sc)
% xlabel('Dimensionless distance (x)')
% ylabel('Concentration (mol/m^3)')
% title('C_{s,n} variation with distance at various times')
% 
% % Electrolyte concentration variation with time
% figure
% plot(t,[C2_p, C2_m,C2_n]*Ce_sc)
% xlabel('Time, sec')
% ylabel('Concentration (mol/m^3)')
% title('Electrolyte concentration with time at various positions')
% 
% % Electrolyte potential variation with time
% figure
% plot(t,[P2_p,P2_m,P2_n]*phie_sc)
% xlabel('Time, sec')
% ylabel('Potential (V)')
% title('Electrolyte potential with time at various positions')



