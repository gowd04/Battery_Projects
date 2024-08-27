function dydt = battery_pde_model_DAE(~,y,Index_chdisch,I_app,Vcell_cc,R_sei, R_plating,Cs_max_p,Cs_max_n,epse_m)

parameters3

C1_p = y(1:N_p);                                     %42
C1_n = y(N_p+1:N_p+N_n);                             %42
C2_p = y(N_p+N_n+1:2*N_p+N_n);                       %42
C2_m = y(2*N_p+N_n+1:2*N_p+N_m+N_n);                 %14
C2_n = y(2*N_p+N_m+N_n+1:2*N_p+N_m+2*N_n);           %42
P1_p = y(2*N_p+N_m+2*N_n+1:3*N_p+N_m+2*N_n);         %42
P1_n = y(3*N_p+N_m+2*N_n+1:3*N_p+N_m+3*N_n);         %42
P2_p = y(3*N_p+N_m+3*N_n+1:4*N_p+N_m+3*N_n);         %42
P2_m = y(4*N_p+N_m+3*N_n+1:4*N_p+2*N_m+3*N_n);       %14
P2_n = y(4*N_p+2*N_m+3*N_n+1:4*N_p+2*N_m+4*N_n);     %42
C3_p = y(4*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+4*N_n);     %42
C3_n = y(5*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+5*N_n);     %42
T_p  = y(5*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+5*N_n);     %42
T_m  = y(6*N_p+2*N_m+5*N_n+1:6*N_p+3*N_m+5*N_n);     %14
T_n  = y(6*N_p+3*N_m+5*N_n+1:6*N_p+3*N_m+6*N_n);     %42

kp_arr  = kp_ref*exp((-Ea_rp/R)*((1./(T_p*Ta))-(1/Tref))); % rate constant for positive electrode
kn1_arr = kn1_ref*exp((-Ea_rn/R)*((1./(T_n*Ta))-(1/Tref)));
kn2_arr = kn2_ref*exp((-Ea_rn/R)*((1./(T_p*Ta))-(1/Tref)));
kn_sei_arr  = kn_sei_ref*exp((-Ea_sei/R)*((1./(T_n*Ta))-(1/Tref)));
kn_plating_arr = kn_plating_ref*exp((-Ea_plating/R)*((1./(T_n*Ta))-(1/Tref)));

kp = mean(kp_arr);
kn1 = mean(kn1_arr);
kn2 = mean(kn2_arr);
kn_sei = mean(kn_sei_arr);
kn_plating = mean(kn_plating_arr);

% kp = kp_ref;
% kn1 = kn1_ref;
% kn2 = kn2_ref;
% kn_sei = kn_sei_ref;
% kn_plating = kn_plating_ref;

if Index_chdisch == 3
    I_app = y(7*N_p+3*N_m+8*N_n+1);
end

Jpi       = y(6*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+6*N_n).*abs(I_app/S_p);
Jni_des   = y(7*N_p+3*N_m+6*N_n+1:7*N_p+3*N_m+7*N_n).*abs(I_app/S_n);

if Index_chdisch ~= 1
    Jni_sei   = y(7*N_p+3*N_m+7*N_n+1:7*N_p+3*N_m+8*N_n)*abs(I_app/S_n);
end

Jpi0          = (kp*F*((C3_p*Csep_sc).^(1-alphaA)).*((C2_p*Ce_sc).^alphaA));
Jni0_des      = (((kn1)^(1-alphaA))* ((kn2)^(1-alphaA)) * F*((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));
Jni0_sei      = (((kn1)^(1-alphaA))* ((kn2 + kn_sei)^(1-alphaA))* F *((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));
% Jni0_sei      = (((kn1).^(1-alphaA))* ((kn2 + kn_sei + kn_plating)^(1-alphaA)) * F *((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));
Jni0_plating  = (((kn1).^(1-alphaA))* ((kn2 + kn_sei + kn_plating)^(1-alphaA)) * F *((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));

if Index_chdisch == 4
    Jni_plating   = y(7*N_p+3*N_m+8*N_n+1:7*N_p+3*N_m+9*N_n)*abs(I_app/S_n);
end

%
%% Diffusivities & Conductivities - Positive electrode
Ds_p = Ds_p_ref*exp((-Ea_dp/R)*((1./(T_p*Ta))-(1/Tref))); % solid phase diffusion coefficient
De_p = De_p_ref*exp((-Ea_e/R)*((1./(T_p*Ta))-(1/Tref)));
% kappae_p = kappae_p*ones(N_p);
kappae_p = kappae_p_ref*exp((-Ea_k/R)*((1./(T_p*Ta))-(1/Tref)));

%% Diffusivities & Conductivities - Separator
De_m = De_m_ref*exp((-Ea_e/R)*((1./(T_m*Ta))-(1/Tref)));
% kappae_m = kappae_m*ones(N_m);
kappae_m = kappae_m_ref*exp((-Ea_k/R)*((1./(T_m*Ta))-(1/Tref)));
 
%% Diffusivities & Conductivities - Negative electrode
Ds_n     = Ds_n_ref*exp((-Ea_dn/R)*((1./(T_n*Ta))-(1/Tref))); % solid phase diffusion coefficient
De_n     = De_n_ref*exp((-Ea_e/R)*((1./(T_n*Ta))-(1/Tref)));
% kappae_n = kappae_n*ones(N_n);
kappae_n = kappae_n_ref*exp((-Ea_k/R)*((1./(T_n*Ta))-(1/Tref)));

%% Diffusivities & Conductivities - Positive electrode
% Ds_p = Ds_p_ref*exp((-Ea_dp/R)*((1./(T_p*Ta))-(1/Ta))); % solid phase diffusion coefficient
% De_p = 1e-4*(10.^((-4.43*(54./((T_p*Ta)-229-(5e-3*(C2_p*Ce_sc))))) - (0.22e-3*(C2_p*Ce_sc))));
% kappae_p = 1e-4*(C2_p*Ce_sc).*((-10.5 + (0.668e-3*(C2_p*Ce_sc)) + (0.494e-6*((C2_p*Ce_sc).^2)) + (0.074*T_p*Ta) -...
%      (1.78e-5*(C2_p*Ce_sc).*(T_p*Ta)) - (8.86e-10*((C2_p*Ce_sc).^2).*(T_p*Ta) - 6.96e-5*((T_p*Ta).^2)) + 2.8e-8*((C2_p*Ce_sc)).*((T_p*Ta).^2)).^2);
% 
% %% Diffusivities & Conductivities - Separator
% 
% De_m = 1e-4*(10.^((-4.43*(54./((T_m*Ta)-229-(5e-3*(C2_m*Ce_sc))))) - (0.22e-3*(C2_m*Ce_sc))));
% kappae_m = 1e-4*(C2_m*Ce_sc).*((-10.5 + (0.668e-3*(C2_m*Ce_sc)) + (0.494e-6*((C2_m*Ce_sc).^2)) + (0.074*T_m*Ta) -...
%     (1.78e-5*(C2_m*Ce_sc).*(T_m*Ta)) - (8.86e-10*((C2_m*Ce_sc).^2).*(T_m*Ta)) - 6.96e-5*((T_m*Ta).^2) + 2.8e-8*((C2_m*Ce_sc)).*((T_m*Ta).^2)).^2);
% 
% %% Diffusivities & Conductivities - Negative electrode
% Ds_n = Ds_n_ref*exp((-Ea_dn/R)*((1./(T_n*Ta))-(1/Ta))); % solid phase diffusion coefficient
% De_n = 1e-4*(10.^((-4.43*(54./((T_n*Ta)-229-(5e-3*(C2_n*Ce_sc))))) - (0.22e-3*(C2_n*Ce_sc))));
% kappae_n = 1e-4*(C2_n*Ce_sc).*((-10.5 + (0.668e-3*(C2_n*Ce_sc)) + (0.494e-6*((C2_n*Ce_sc).^2)) + (0.074*T_n*Ta) -...
%     (1.78e-5*(C2_n*Ce_sc).*(T_n*Ta)) - (8.86e-10*((C2_n*Ce_sc).^2).*(T_n*Ta) - 6.96e-5*((T_n*Ta).^2)) + 2.8e-8*((C2_n*Ce_sc)).*((T_n*Ta).^2)).^2);

%% Current density, theta and OCV
Jptot = I_app*ones(N_p,1)/(as_p*S_p*L_p);
Jntot = -I_app*ones(N_n,1)/(as_n*S_n*L_n);

theta_p    = C3_p*Csep_sc/Cs_max_p;
theta_n    = C3_n*Csen_sc/Cs_max_n;

% Equation obtained from 'Modeling of lithium plating and lithium stripping in lithium-ion batteries'

Uref_p = 6.0826 - 6.9922*theta_p + 7.1062*(theta_p.^2) - 0.54549e-4*exp((124.23 * theta_p) - 114.2593) - (2.5947 * theta_p.^3);

Uref_n = 0.6379 + 0.5416*exp(-305.5309*theta_n) + 0.044*tanh(-(theta_n - 0.1958)/0.1088) - 0.1978*tanh((theta_n - 1.0571)/0.0854) - 0.6875*tanh((theta_n + 0.0117)/0.0529)...
    - 0.0175*tanh((theta_n - 0.5692)/0.0875);

% Equation obtained by fitting curves from 'Analysis of Capacity Fade from Entropic Heat
% Coefficient of Li[NixCoyMnz]O2/Graphite Lithium Ion Battery'

Ent_coeff_p = -4.8833e02 *theta_p.^15 + 3.7717e03 *theta_p.^14 - 1.3181e04*theta_p.^13 + 2.7550e04*theta_p.^12 ...
        -3.8330e04*theta_p.^11 + 3.7379e04*theta_p.^10 - 2.6196e04*theta_p.^9 + 1.3304e04*theta_p.^8 - 4.8754e03*theta_p.^7 ...
        + 1.2678e03*theta_p.^6 - 2.2710e02*theta_p.^5 + 2.6738e01*theta_p.^4 - 1.9261*theta_p.^3 ...
        + 7.5184e-02*theta_p.^2 - 1.3270e-03*theta_p + 7.5369e-05;

Ent_coeff_n = -1.9306e03*theta_n.^15 + 1.4344e04*theta_n.^14 - 4.7972e04*theta_n.^13 + 9.5423e04*theta_n.^12 ...
        -1.2568e05*theta_n.^11 + 1.1546e05*theta_n.^10 - 7.5903e04*theta_n.^9 + 3.6017e04*theta_n.^8 - 1.2283e04*theta_n.^7 ...
        + 2.9583e03*theta_n.^6 - 4.8684e02*theta_n.^5 + 5.1903e01*theta_n.^4 - 3.3039*theta_n.^3 ...
        + 1.0959e-01*theta_n.^2 + 5.6345e-04*theta_n - 3.6876e-04;

Ui_p = Uref_p + ((T_p - 1).*(Ent_coeff_p));
Ui_n = Uref_n + ((T_n - 1).*(Ent_coeff_n));

etap = P1_p*phisp_sc - P2_p*phie_sc - Ui_p; %Overpotential

if Index_chdisch == 1
    etan = P1_n*phisn_sc - P2_n*phie_sc - Ui_n - Jni_des.*(R_sei' + R_plating');
elseif Index_chdisch == 4
    etan = P1_n*phisn_sc - P2_n*phie_sc - Ui_n - Jni_plating.*(R_sei' + R_plating');
else
    etan = P1_n*phisn_sc - P2_n*phie_sc - Ui_n - Jni_sei.*(R_sei' + R_plating');
end

%% 1. Derivatives of solid phase concentration
dC1p_dt = -(as_p*Jptot*t_sc)/(F*epss_p*Csp_sc);
if Index_chdisch == 1
    dC1n_dt = -(as_n*Jntot*t_sc)/(F*epss_n*Csn_sc);
elseif Index_chdisch == 4
    dC1n_dt = -(as_n*Jntot*t_sc.*Jni0_des)./(F*epss_n*Csn_sc*Jni0_plating);
% else
%     dC1n_dt = -(as_n*Jntot*t_sc.*Jni0_des)./(F*epss_n*Csn_sc*Jni0_plating);
else
    dC1n_dt = -(as_n*Jntot*t_sc.*Jni0_des)./(F*epss_n*Csn_sc*Jni0_sei);
end

%% 2. Derivatives of electrolyte phase conc in positive electrode
dC2p_dt(1) = C2_p(2) - C2_p(1);

for kdx = 2:N_p-1
    dC2p_dt(kdx,1) = ((De_p(kdx)*(epse_p^(Brugg_p-1))*t_sc/((x_sc*Delx_p)^2))*(C2_p(kdx+1) - 2*C2_p(kdx) + C2_p(kdx-1)))...
        +(((1-transf_p(kdx))*t_sc*as_p*Jpi(kdx))/(F*epse_p*Ce_sc));
end

dC2p_dt(N_p) = (((De_m(1)*(epse_m^(Brugg_m))*Delx_p)/(De_p(N_p)*(epse_p^(Brugg_p))*Delx_m))*(C2_m(1) - C2_p(N_p)))  - (C2_p(N_p) - C2_p (N_p - 1));

%% 3. Derivatives of electrolyte phase conc in membrane separator

dC2m_dt(1) = ((De_m(1)*(epse_m^(Brugg_m-1))*t_sc/((x_sc*Delx_m)^2))*(C2_m(2) - 2*C2_m(1) + C2_p(N_p)));
for kdx = 2:N_m - 1
    dC2m_dt(kdx,1) = ((De_m(kdx)*(epse_m^(Brugg_m-1))*t_sc/((x_sc*Delx_m)^2))*(C2_m(kdx+1) - 2*C2_m(kdx) + C2_m(kdx-1)));
end
dC2m_dt(N_m) = ((De_m(N_m)*(epse_m^(Brugg_m-1))*t_sc/((x_sc*Delx_m)^2))*(C2_n(1) - 2*C2_m(N_m) + C2_m(N_m-1)));

%% 4. Derivatives of electrolyte phase conc in negative electrode
dC2n_dt(1) = (((De_m(1)*(epse_m^(Brugg_m))*Delx_n)/(De_n(1)*(epse_n^(Brugg_n))*Delx_m))*(C2_n(1) - C2_m(N_m)))  - (C2_n(2) - C2_n (1));

if Index_chdisch ==1
    for kdx = 2:N_n-1
        dC2n_dt(kdx,1) = ((De_n(kdx)*(epse_n^(Brugg_n-1))*t_sc/((x_sc*Delx_n)^2))*(C2_n(kdx+1) - 2*C2_n(kdx) + C2_n(kdx-1)))...
            +(((1-transf_n(kdx))*t_sc*as_n*Jni_des(kdx))/(F*epse_n*Ce_sc));
    end
elseif Index_chdisch == 4
    for kdx = 2:N_n-1
        dC2n_dt(kdx,1) = ((De_n(kdx)*(epse_n^(Brugg_n-1))*t_sc/((x_sc*Delx_n)^2))*(C2_n(kdx+1) - 2*C2_n(kdx) + C2_n(kdx-1)))...
            +(((1-transf_n(kdx))*t_sc*as_n*Jni_plating(kdx))/(F*epse_n*Ce_sc));
    end
else
    for kdx = 2:N_n-1
        dC2n_dt(kdx,1) = ((De_n(kdx)*(epse_n^(Brugg_n-1))*t_sc/((x_sc*Delx_n)^2))*(C2_n(kdx+1) - 2*C2_n(kdx) + C2_n(kdx-1)))...
            +(((1-transf_n(kdx))*t_sc*as_n*Jni_sei(kdx))/(F*epse_n*Ce_sc));
    end
end
dC2n_dt(N_n) = C2_n(N_n) - C2_n(N_n - 1);

%% 5. Derivatives of solid phase potential in positive electrode

dP1p_dt(1) = (P1_p(2) - P1_p(1)) +((Delx_p)*I_app*x_sc/(S_p*phisp_sc*sigmas_p*(epss_p^(Brugg_p))));

for kdx = 2:N_p-1
    dP1p_dt(kdx,1) = -((Delx_p^2)*as_p*Jpi(kdx)*(x_sc^2)/(phisp_sc*sigmas_p*(epss_p^(Brugg_p)))) + (P1_p(kdx+1) - 2*P1_p(kdx) + P1_p(kdx-1)) ;
end
dP1p_dt(N_p) = P1_p(N_p) - P1_p(N_p - 1);

%% 6. Derivatives of solid phase potential in negative electrode

dP1n_dt(1) = P1_n(2) - P1_n(1);

for kdx = 2:N_n-1
    dP1n_dt(kdx,1) = -((Delx_n^2)*as_n*Jni_des(kdx)*(x_sc^2)/(phisn_sc*sigmas_n*(epss_n^(Brugg_n)))) + ((P1_n(kdx+1) - 2*P1_n(kdx) + P1_n(kdx-1)));
end

dP1n_dt(N_n) = (P1_n(N_n) - P1_n(N_n -1)) + ((Delx_n)*I_app*x_sc/(S_n*phisn_sc*sigmas_n*(epss_n^(Brugg_n))))  ;

%% 7. Derivatives of electrolyte phase potential in positive electrode

dP2p_dt(1) = P2_p(2) - P2_p(1);

for kdx = 2:N_p-1
    dP2p_dt(kdx,1) = ((P2_p(kdx+1) - 2*P2_p(kdx) + P2_p(kdx-1))) + ((2*R*T_p(kdx)*Ta*(transf_p(kdx) - 1)*(log(C2_p(kdx+1)) - 2*log(C2_p(kdx)) +log(C2_p(kdx-1))))/(phie_sc*F)) ...
        + (((Delx_p^2)*as_p*Jpi(kdx)*(x_sc^2))/(phie_sc*kappae_p(kdx)*(epse_p^(Brugg_p))));
end

dP2p_dt(N_p) =  ((kappae_m(1)*(epse_m^(Brugg_m))*Delx_p)/(kappae_p(N_p)*(epse_p^(Brugg_p))*Delx_m))*(P2_m(1) - P2_p(N_p)) - (P2_p(N_p) - P2_p(N_p - 1));

%% 8. Derivatives of electrolyte phase potential in membrane separator

dP2m_dt(1) = (P2_m(2) - 2*P2_m(1) + P2_p(N_n)) + ((2*R*T_m(1)*Ta*(transf_m(1) - 1)*(log(C2_m(2)) - 2*log(C2_m(1)) +log(C2_p(N_p))))/(phie_sc*F)) ;

for kdx = 2:N_m-1
    dP2m_dt(kdx,1) = (P2_m(kdx+1) - 2*P2_m(kdx) + P2_m(kdx-1)) + ((2*R*T_m(kdx)*Ta*(transf_m(kdx) - 1)*(log(C2_m(kdx+1)) - 2*log(C2_m(kdx)) +log(C2_m(kdx-1))))/(phie_sc*F));
end

dP2m_dt(N_m) = (P2_n(1) - 2*P2_m(N_m) + P2_m(N_m-1)) + ((2*R*T_m(N_m)*Ta*(transf_m(N_m) - 1)*(log(C2_n(1)) - 2*log(C2_m(N_m)) +log(C2_m(N_m-1))))/(phie_sc*F));

%% 9. Derivatives of electrolyte phase potential in negative electrode

dP2n_dt(1) =  ((kappae_n(1)*(epse_m^(Brugg_m))*Delx_n)/(kappae_n(1)*(epse_n^(Brugg_n))*Delx_m))*(P2_n(1) - P2_m(N_m)) - (P2_n(2) - P2_n(1));

if Index_chdisch == 1
    for kdx = 2:N_n-1
        dP2n_dt(kdx,1) = (P2_n(kdx+1) - 2*P2_n(kdx) + P2_n(kdx-1)) + ((2*R*T_n(kdx)*Ta*(transf_n(kdx) - 1)*(log(C2_n(kdx+1)) - 2*log(C2_n(kdx)) +log(C2_n(kdx-1))))/(phie_sc*F)) ...
            + (((Delx_n^2)*as_n*Jni_des(kdx)*(x_sc^2))/(phie_sc*kappae_n(kdx)*(epse_n^(Brugg_n))));
    end
elseif Index_chdisch == 4
    for kdx = 2:N_n-1
        dP2n_dt(kdx,1) = (P2_n(kdx+1) - 2*P2_n(kdx) + P2_n(kdx-1)) + ((2*R*T_n(kdx)*Ta*(transf_n(kdx) - 1)*(log(C2_n(kdx+1)) - 2*log(C2_n(kdx)) +log(C2_n(kdx-1))))/(phie_sc*F)) ...
            + (((Delx_n^2)*as_n*Jni_plating(kdx)*(x_sc^2))/(phie_sc*kappae_n(kdx)*(epse_n^(Brugg_n))));
    end
else
    for kdx = 2:N_n-1
        dP2n_dt(kdx,1) = (P2_n(kdx+1) - 2*P2_n(kdx) + P2_n(kdx-1)) + ((2*R*T_n(kdx)*Ta*(transf_n(kdx) - 1)*(log(C2_n(kdx+1)) - 2*log(C2_n(kdx)) +log(C2_n(kdx-1))))/(phie_sc*F)) ...
            + (((Delx_n^2)*as_n*Jni_sei(kdx)*(x_sc^2))/(phie_sc*kappae_n(kdx)*(epse_n^(Brugg_n))));
    end
end
dP2n_dt(N_n) = P2_n(N_n);

%% 10. Derivatives of concentrations at solid electrolyte interphase

dC3p_dt = C3_p - C1_p*Csp_sc/Csep_sc +(Jptot*Rs_p./(5*F*Ds_p*Csep_sc));

if Index_chdisch == 1
    dC3n_dt = C3_n - C1_n*Csn_sc/Csen_sc +(Jntot.*Rs_n./(5*F*Ds_n.*Csen_sc));
elseif Index_chdisch == 4
    dC3n_dt = C3_n - C1_n*Csn_sc/Csen_sc +(Jntot.*Rs_n.*Jni0_des./(5*F*Ds_n.*Csen_sc.*Jni0_plating));
% else
%     dC3n_dt = C3_n - C1_n*Csn_sc/Csen_sc +(Jntot.*Rs_n.*Jni0_des./(5*F*Ds_n.*Csen_sc.*Jni0_plating));
else
    dC3n_dt = C3_n - C1_n*Csn_sc/Csen_sc +(Jntot.*Rs_n.*Jni0_des./(5*F*Ds_n.*Csen_sc.*Jni0_sei));
end

%% 11. Energy balance
    dTdt = zeros(N_p+N_m+N_n,1);
    for kdx = 1:(N_p+N_m+N_n)
        if kdx <= N_p && kdx >= 2 %% positive electrode
            if kdx < N_p
    %             Ui_p(kdx) = Uref_p(kdx) + ((T_p(kdx)-Ta)*Ent_coeff_p(kdx));
                Qrxn = as_p*Jpi(kdx)*((P1_p(kdx)*phisp_sc) - (P2_p(kdx)*phie_sc) - Ui_p(kdx));
                Qrev = as_p*Jpi(kdx)*T_p(kdx)*Ent_coeff_p(kdx);
                Qohm = (sigmas_p*(epss_p^(Brugg_p))*(((P1_p(kdx +1) - P1_p(kdx))/Delx_p)^2)*((phisp_sc/x_sc)^2)) + (kappae_p(kdx)*(epse_p^(Brugg_p))*(((P2_p(kdx +1) - P2_p(kdx))/Delx_p)^2)*((phie_sc/x_sc)^2)) + ...
                    ((2*R*kappae_p(kdx)*(epse_p^(Brugg_p))*T_p(kdx)*Ta*(1-transf_p(kdx))*(((C2_p(kdx+1))-(C2_p(kdx)))/Delx_p)*((P2_p(kdx +1) - P2_p(kdx))/Delx_p)*phie_sc)/(F*C2_p(kdx)*(x_sc^2)));
                dTdt(kdx,1) = ((lambda_p*t_sc)/(rho_p*Cp_p*(x_sc^2)))*((T_p(kdx+1) - 2*T_p(kdx) + T_p(kdx-1))/(Delx_p^2)) + (t_sc*Qrxn)/(rho_p*Cp_p*Ta) + (t_sc*Qrev)/(rho_p*Cp_p*Ta) + (t_sc*Qohm)/(rho_p*Cp_p*Ta);
                
            else
    %             Ui_p(kdx) = Uref_p(kdx) + ((T_p(kdx)-Ta)*Ent_coeff_p(kdx));
                Qrxn = as_p*Jpi(kdx)*((P1_p(kdx)*phisp_sc) - (P2_p(kdx)*phie_sc) - Ui_p(kdx));
                Qrev = as_p*Jpi(kdx)*T_p(kdx)*Ent_coeff_p(kdx);
                Qohm = (sigmas_p*(epss_p^(Brugg_p))*(((P1_p(kdx) - P1_p(kdx-1))/Delx_p)^2)*((phisp_sc/x_sc)^2)) + (kappae_p(kdx)*(epse_p^(Brugg_p))*(((P2_m(1) - P2_p(kdx))/Delx_p)^2)*((phie_sc/x_sc)^2)) + ...
                    ((2*R*kappae_p(kdx)*(epse_p^(Brugg_p))*T_p(kdx)*Ta*(1-transf_p(kdx))*((C2_m(1)-C2_p(kdx))/Delx_p)*((P2_m(1) - P2_p(kdx))/Delx_p)*phie_sc)/(F*C2_p(kdx)*(x_sc^2)));
                dTdt(kdx,1) = ((lambda_p*t_sc)/(rho_p*Cp_p*(x_sc^2)))*((T_m(1) - 2*T_p(kdx) + T_p(kdx-1))/(Delx_p^2)) + (t_sc*Qrxn)/(rho_p*Cp_p*Ta) + (t_sc*Qrev)/(rho_p*Cp_p*Ta) + (t_sc*Qohm)/(rho_p*Cp_p*Ta);
                
            end
        end
        
        if kdx > N_p && kdx <= N_p+N_m  %%separator
            if kdx == N_p+1 %interphase at positive electrode
                Qohm = (kappae_m(kdx-N_p)*(epse_m^(Brugg_m))*(((P2_m(kdx+1 - N_p) - P2_m(kdx - N_p))/Delx_m)^2)*((phie_sc/x_sc)^2)) +...
                    ((2*R*kappae_m(kdx-N_p)*(epse_m^(Brugg_m))*T_m(kdx - N_p)*Ta*(1-transf_m(kdx-N_p))*(((C2_m(kdx+1 - N_p)) - (C2_m(kdx - N_p)))/Delx_m)*((P2_m(kdx +1-N_p) - P2_m(kdx - N_p))/Delx_m)*phie_sc)/(F*C2_m(kdx-N_p)*(x_sc^2)));
                dTdt(kdx,1) = ((lambda_m*t_sc)/(rho_m*Cp_m*(x_sc^2)))*(((T_m(kdx+1-N_p)) - 2*(T_m(kdx - N_p)) + T_p(N_p))/(Delx_m^2)) + (t_sc*Qohm)/(rho_m*Cp_m*Ta);
                
            elseif kdx == N_p+N_m
                Qohm = (kappae_m(kdx-N_p)*(epse_m^(Brugg_m))*(((P2_n(1) - P2_m(kdx - N_p))/Delx_m)^2)*((phie_sc/x_sc)^2)) +...
                    ((2*R*kappae_m(kdx-N_p)*(epse_m^(Brugg_m))*T_m(kdx - N_p)*Ta*(1-transf_m(kdx-N_p))*((C2_n(1) - (C2_m(kdx - N_p)))/Delx_m)*((P2_n(1) - P2_m(kdx - N_p))/Delx_m)*phie_sc)/(F*C2_m(kdx-N_p)*(x_sc^2)));
                dTdt(kdx,1) = ((lambda_m*t_sc)/(rho_m*Cp_m*(x_sc^2)))*((T_n(1) - 2*(T_m(kdx - N_p)) + T_m(kdx-1-N_p))/(Delx_m^2)) + ((t_sc*Qohm)/(rho_m*Cp_m*Ta));
                
            else
                Qohm = (kappae_m(kdx-N_p)*(epse_m^(Brugg_m))*(((P2_m(kdx+1 - N_p) - P2_m(kdx - N_p))/Delx_m)^2)*((phie_sc/x_sc)^2)) +...
                    ((2*R*kappae_m(kdx-N_p)*(epse_m^(Brugg_m))*T_m(kdx - N_p)*Ta*(1-transf_m(kdx-N_p))*(((C2_m(kdx+1 - N_p)) - (C2_m(kdx - N_p)))/Delx_m)*((P2_m(kdx +1-N_p) - P2_m(kdx-N_p))/Delx_m)*phie_sc)/(F*C2_m(kdx-N_p)*(x_sc^2)));
                dTdt(kdx,1) = ((lambda_m*t_sc)/(rho_m*Cp_m*(x_sc^2)))*((T_m(kdx+1-N_p) - 2*(T_m(kdx - N_p)) + T_m(kdx-1-N_p))/(Delx_m^2)) + (t_sc*Qohm)/(rho_m*Cp_m*Ta);
            end
        end
        
        if kdx > N_p+N_m && kdx < N_n+N_m+N_p %% negative electrode
            if kdx == N_p+N_m+1
                if Index_chdisch == 1
                    Jn = Jni_des;
                elseif Index_chdisch == 4
                    Jn = Jni_plating;
                else
                    Jn = Jni_sei;
                end
    %             Ui_n(kdx-(N_p+N_m)) = Uref_n(kdx-(N_p+N_m)) + ((T_n(kdx-(N_p+N_m))-Ta)*Ent_coeff_n(kdx-(N_p+N_m)));
                Qrxn = as_n*Jn(kdx-(N_p+N_m))*((P1_n(kdx-(N_p+N_m))*phisn_sc) - (P2_n(kdx-(N_p+N_m))*phie_sc) - Ui_n(kdx-(N_p+N_m)));
                Qrev = as_n*Jn(kdx-(N_p+N_m))*T_n(kdx-(N_p+N_m))*Ent_coeff_n(kdx-(N_p+N_m));
                Qohm = ((sigmas_n*(epss_n^(Brugg_n))*(((P1_n(kdx +1-(N_p+N_m)) - P1_n(kdx-(N_p+N_m)))/Delx_n)^2))*((phisn_sc/x_sc)^2)) + ((kappae_n(kdx-(N_p+N_m))*(epse_n^(Brugg_n))*(((P2_n(kdx +1-(N_p+N_m)) - P2_n(kdx-(N_p+N_m)))/Delx_n)^2))*((phie_sc/x_sc)^2)) + ...
                    ((2*R*kappae_n(kdx-(N_p+N_m))*(epse_n^(Brugg_n))*T_n(kdx-(N_p+N_m))*Ta*(1-transf_n(kdx-(N_p+N_m)))*(((C2_n(kdx+1-(N_p+N_m)))-(C2_n(kdx-(N_p+N_m))))/Delx_n)*((P2_n(kdx +1-(N_p+N_m)) - P2_n(kdx-(N_p+N_m)))/Delx_n)*phie_sc)/(F*C2_n(kdx-(N_p+N_m))*(x_sc^2)));
                dTdt(kdx,1) = ((lambda_n*t_sc)/(rho_n*Cp_n*(x_sc^2)))*((T_n(kdx+1-(N_p+N_m)) - 2*T_n(kdx-(N_p+N_m)) + T_m(N_m))/(Delx_n^2)) + (t_sc*Qrxn)/(rho_n*Cp_n*Ta) + (t_sc*Qrev)/(rho_n*Cp_n*Ta) + (t_sc*Qohm)/(rho_n*Cp_n*Ta);
            else
                if Index_chdisch == 1
                    Jn = Jni_des;
                elseif Index_chdisch == 4
                    Jn = Jni_plating;
                else
                    Jn = Jni_sei;
                end
    %             Ui_n(kdx-(N_p+N_m)) = Uref_n(kdx-(N_p+N_m)) + ((T_n(kdx-(N_p+N_m))-Ta)*Ent_coeff_n(kdx-(N_p+N_m)));
                Qrxn = as_n*Jn(kdx-(N_p+N_m))*(P1_n(kdx-(N_p+N_m))*phisn_sc - P2_n(kdx-(N_p+N_m))*phie_sc - Ui_n(kdx-(N_p+N_m)));
                Qrev = as_n*Jn(kdx-(N_p+N_m))*T_n(kdx-(N_p+N_m))*Ent_coeff_n(kdx-(N_p+N_m));
                Qohm = (sigmas_n*(epss_n^(Brugg_n))*(((P1_n(kdx +1-(N_p+N_m)) - P1_n(kdx-(N_p+N_m)))/Delx_n)^2)*((phisn_sc/x_sc)^2)) + (kappae_n(kdx-(N_p+N_m))*(epse_n^(Brugg_n))*(((P2_n(kdx +1-(N_p+N_m)) - P2_n(kdx-(N_p+N_m)))/Delx_n)^2)*((phie_sc/x_sc)^2)) + ...
                    ((2*R*kappae_n(kdx-(N_p+N_m))*(epse_n^(Brugg_n))*T_n(kdx-(N_p+N_m))*Ta*(1-transf_n(kdx-(N_p+N_m)))*(((C2_n(kdx+1-(N_p+N_m)))-(C2_n(kdx-(N_p+N_m))))/Delx_n)*((P2_n(kdx +1-(N_p+N_m)) - P2_n(kdx-(N_p+N_m)))/Delx_n)*phie_sc)/(F*C2_n(kdx-(N_p+N_m))*(x_sc^2)));
                
                dTdt(kdx,1) = ((lambda_n*t_sc)/(rho_n*Cp_n*(x_sc^2)))*((T_n(kdx+1-(N_p+N_m)) - 2*T_n(kdx-(N_p+N_m)) + T_n(kdx-1-(N_p+N_m)))/(Delx_n^2)) + (t_sc*Qrxn)/(rho_n*Cp_n*Ta)+ (t_sc*Qrev)/(rho_n*Cp_n*Ta)  + (t_sc*Qohm)/(rho_n*Cp_n*Ta);
            end
        end
        
        if kdx == 1
            dTdt(kdx,1) = ((-T_p(kdx) + T_p(kdx+1))/(Delx_p*x_sc) + ((h*(1-T_p(kdx)))/lambda_p));
        end
        if kdx == N_n+N_m+N_p
            dTdt(kdx,1) = ((T_n(end) - T_n(end-1))/(Delx_n*x_sc) + ((h*(T_n(end)-1))/lambda_n));
        end
    end

%% 12. Current density equations

    dJpidt = (Jpi - Jpi0.*(exp((alphaA*F*etap)./(R*T_p*Ta)) - exp(-((1-alphaA)*F*etap)./(R*T_p*Ta))))/abs(I_app/S_p);
    
    dJnidt_des = (Jni_des - Jni0_des.*(exp((alphaA*F*etan)./(R*T_n*Ta)) - exp(-((1-alphaA)*F*etan)./(R*T_n*Ta))))/abs(I_app/S_p);
    
    if Index_chdisch ~= 1
        dJnidt_sei = (Jni_sei - Jni0_sei.*(exp((alphaA*F*etan)./(R*T_n*Ta)) - exp(-((1-alphaA)*F*etan)./(R*T_n*Ta))))/abs(I_app/S_p);
    end
    if Index_chdisch == 4
        dJnidt_plating = (Jni_plating - Jni0_plating.*(exp((alphaA*F*etan)./(R*T_n*Ta)) - exp(-((1-alphaA)*F*etan)./(R*T_n*Ta))))/abs(I_app/S_p);
    end

    if Index_chdisch == 1
        
        %% Final expression
        dydt = [dC1p_dt; dC1n_dt; dC2p_dt; dC2m_dt; dC2n_dt; ...
            dP1p_dt; dP1n_dt; dP2p_dt; dP2m_dt; dP2n_dt;dC3p_dt;dC3n_dt;dTdt;dJpidt;dJnidt_des];
        
    elseif Index_chdisch == 3
        
        %% voltage of cell
        dIapp_dt = P1_p(1)* phisp_sc - P1_n(end)*phisn_sc - Vcell_cc;
        
        %% Final expression
        dydt = [dC1p_dt; dC1n_dt; dC2p_dt; dC2m_dt; dC2n_dt; ...
            dP1p_dt; dP1n_dt; dP2p_dt; dP2m_dt; dP2n_dt; dC3p_dt; dC3n_dt;dTdt; dJpidt; dJnidt_des; dJnidt_sei; dIapp_dt];
    else
        
        %% Final expression
        dydt = [dC1p_dt; dC1n_dt; dC2p_dt; dC2m_dt; dC2n_dt; ...
            dP1p_dt; dP1n_dt; dP2p_dt; dP2m_dt; dP2n_dt; dC3p_dt; dC3n_dt;dTdt; dJpidt; dJnidt_des; dJnidt_sei; dJnidt_plating];
        
    end

end

