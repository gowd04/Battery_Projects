function dydt = battery_pde_model_DAE(~,y,Index_chdisch,I_app,Vcell_cc,R_sei, R_plating,Cs_max_p,Cs_max_n,epse_m)
    
    parameters2
    
    C1_p = y(1:N_p);
    C1_n = y(N_p+1:N_p+N_n);
    C2_p = y(N_p+N_n+1:2*N_p+N_n);
    C2_m = y(2*N_p+N_n+1:2*N_p+N_m+N_n);
    C2_n = y(2*N_p+N_m+N_n+1:2*N_p+N_m+2*N_n);
    P1_p = y(2*N_p+N_m+2*N_n+1:3*N_p+N_m+2*N_n);
    P1_n = y(3*N_p+N_m+2*N_n+1:3*N_p+N_m+3*N_n);
    P2_p = y(3*N_p+N_m+3*N_n+1:4*N_p+N_m+3*N_n);
    P2_m = y(4*N_p+N_m+3*N_n+1:4*N_p+2*N_m+3*N_n);
    P2_n = y(4*N_p+2*N_m+3*N_n+1:4*N_p+2*N_m+4*N_n);
    C3_p = y(4*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+4*N_n);
    C3_n = y(5*N_p+2*N_m+4*N_n+1:5*N_p+2*N_m+5*N_n);
    
    if Index_chdisch == 3
        I_app = y(6*N_p+2*N_m+7*N_n+1);
    end
    Jpi       = y(5*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+5*N_n)*abs(I_app/S_p);
    Jni_des   = y(6*N_p+2*N_m+5*N_n+1:6*N_p+2*N_m+6*N_n)*abs(I_app/S_n);
    
    if Index_chdisch ~= 1
       Jni_sei   = y(6*N_p+2*N_m+6*N_n+1:6*N_p+2*N_m+7*N_n)*abs(I_app/S_n);
    end
    
    Jpi0          = (kp * F*((C3_p*Csep_sc).^(1-alphaA)).*((C2_p*Ce_sc).^alphaA));
    Jni0_des      = (((kn1)^(1-alphaA))* ((kn2)^(1-alphaA)) * F*((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));
    Jni0_sei      = (((kn1)^(1-alphaA))* ((kn2 + kn_sei)^(1-alphaA)) * F*((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));
    
    if Index_chdisch == 4
        Jni0_plating  = (((kn1)^(1-alphaA))* ((kn2 + kn_sei + kn_plating)^(1-alphaA)) * F*((C3_n*Csen_sc).^(1-alphaA)).*((C2_n*Ce_sc).^alphaA));
        Jni_plating   = y(6*N_p+2*N_m+7*N_n+1:6*N_p+2*N_m+8*N_n)*abs(I_app/S_n);        
    end
    
    Jptot = I_app*ones(N_p,1)/(as_p*S_p*L_p);
    Jntot = -I_app*ones(N_n,1)/(as_n*S_n*L_n);    
    
    theta_p    = C3_p*Csep_sc/Cs_max_p;
    theta_n    = C3_n*Csen_sc/Cs_max_n;
    
    Uref_p = 1.654107793108376e+06 *theta_p.^10 -1.249511527828579e+07 *theta_p.^9 + 4.215812681219930e+07*theta_p.^8 -8.365902527257702e+07*theta_p.^7 + ...
        1.081256432503016e+08*theta_p.^6 -9.510080800780240e+07*theta_p.^5 + 5.764438732509381e+07*theta_p.^4 -2.377606148390258e+07*theta_p.^3 + ...
        6.386329451323811e+06*theta_p.^2 -1.008737242075382e+06*theta_p +  7.115616223683314e+04;
            
    Uref_n = (1.19970 + (118.1911 * (theta_n.^0.5)) - (706.0711 * (theta_n)) + (2217.6479 * (theta_n.^1.5)) - (1675.1321 * (theta_n.^2)))./...
        (1.0 + (131.7572 * (theta_n.^0.5)) - (32.1402 * (theta_n)) - (746.8463 * (theta_n.^1.5)) + (15502.9505 * (theta_n.^2)) - (14213.0747 * (theta_n.^2.5)));
    
    etap = P1_p*phisp_sc - P2_p*phie_sc - Uref_p;
    
    if Index_chdisch == 1
        etan = P1_n*phisn_sc - P2_n*phie_sc - Uref_n - Jni_des.*(R_sei' + R_plating');
    elseif Index_chdisch == 4
        etan = P1_n*phisn_sc - P2_n*phie_sc - Uref_n - Jni_plating.*(R_sei' + R_plating');
    else
        etan = P1_n*phisn_sc - P2_n*phie_sc - Uref_n - Jni_sei.*(R_sei' + R_plating');
    end
    
    %% Derivatives of solid phase concentration
    dC1p_dt = -(as_p*Jptot*t_sc)/(F*epss_p*Csp_sc);
    
    if Index_chdisch == 1
        dC1n_dt = -(as_n*Jntot*t_sc)/(F*epss_n*Csn_sc);
    elseif Index_chdisch == 4
        dC1n_dt = -(as_n*Jntot*t_sc.*Jni0_des)./(F*epss_n*Csn_sc*Jni0_plating);
    else
        dC1n_dt = -(as_n*Jntot*t_sc.*Jni0_des)./(F*epss_n*Csn_sc*Jni0_sei);
    end
    
    %% Derivatives of electrolyte phase conc in positive electrode
    dC2p_dt(1) = C2_p(2) - C2_p(1);
    
    for kdx = 2:N_p-1
        dC2p_dt(kdx,1) = ((De_p*(epse_p^(Brugg_p-1))*t_sc/((x_sc*Delx_p)^2))*(C2_p(kdx+1) - 2*C2_p(kdx) + C2_p(kdx-1)))...
            +(((1-transf)*t_sc*as_p*Jpi(kdx))/(F*epse_p*Ce_sc));
    end
    
    dC2p_dt(N_p) = (((De_m*(epse_m^(Brugg_m))*Delx_p)/(De_p*(epse_p^(Brugg_p))*Delx_m))*(C2_m(1) - C2_p(N_p)))  - (C2_p(N_p) - C2_p (N_p - 1));
    
    %% Derivatives of electrolyte phase conc in membrane separator
    
    dC2m_dt(1) = ((De_m*(epse_m^(Brugg_m-1))*t_sc/((x_sc*Delx_m)^2))*(C2_m(2) - 2*C2_m(1) + C2_p(N_p)));
    for kdx = 2:N_m - 1
        dC2m_dt(kdx,1) = ((De_m*(epse_m^(Brugg_m-1))*t_sc/((x_sc*Delx_m)^2))*(C2_m(kdx+1) - 2*C2_m(kdx) + C2_m(kdx-1)));
    end
    dC2m_dt(N_m) = ((De_m*(epse_m^(Brugg_m-1))*t_sc/((x_sc*Delx_m)^2))*(C2_n(1) - 2*C2_m(N_m) + C2_m(N_m-1)));
    
    %% Derivatives of electrolyte phase conc in negative electrode
    dC2n_dt(1) = (((De_m*(epse_m^(Brugg_m))*Delx_n)/(De_n*(epse_n^(Brugg_n))*Delx_m))*(C2_n(1) - C2_m(N_m)))  - (C2_n(2) - C2_n (1));
    
    if Index_chdisch ==1
        for kdx = 2:N_n-1
            dC2n_dt(kdx,1) = ((De_n*(epse_n^(Brugg_n-1))*t_sc/((x_sc*Delx_n)^2))*(C2_n(kdx+1) - 2*C2_n(kdx) + C2_n(kdx-1)))...
                +(((1-transf)*t_sc*as_n*Jni_des(kdx))/(F*epse_n*Ce_sc));
        end   
    elseif Index_chdisch == 4
        for kdx = 2:N_n-1
            dC2n_dt(kdx,1) = ((De_n*(epse_n^(Brugg_n-1))*t_sc/((x_sc*Delx_n)^2))*(C2_n(kdx+1) - 2*C2_n(kdx) + C2_n(kdx-1)))...
                +(((1-transf)*t_sc*as_n*Jni_plating(kdx))/(F*epse_n*Ce_sc));
        end
    else
        for kdx = 2:N_n-1
            dC2n_dt(kdx,1) = ((De_n*(epse_n^(Brugg_n-1))*t_sc/((x_sc*Delx_n)^2))*(C2_n(kdx+1) - 2*C2_n(kdx) + C2_n(kdx-1)))...
                +(((1-transf)*t_sc*as_n*Jni_sei(kdx))/(F*epse_n*Ce_sc));
        end        
    end
    dC2n_dt(N_n) = C2_n(N_n) - C2_n(N_n - 1);
    
    %% Derivatives of solid phase potential in positive electrode
    
    dP1p_dt(1) = (P1_p(2) - P1_p(1)) +((Delx_p)*I_app*x_sc/(S_p*phisp_sc*sigmas_p*(epss_p^(Brugg_p))));
    
    for kdx = 2:N_p-1
        dP1p_dt(kdx,1) = -((Delx_p^2)*as_p*Jpi(kdx)*(x_sc^2)/(phisp_sc*sigmas_p*(epss_p^(Brugg_p)))) + (P1_p(kdx+1) - 2*P1_p(kdx) + P1_p(kdx-1)) ;
    end
    dP1p_dt(N_p) = P1_p(N_p) - P1_p(N_p - 1);
    
    %% Derivatives of solid phase potential in negative electrode
    
    dP1n_dt(1) = P1_n(2) - P1_n(1);
    
    for kdx = 2:N_n-1
        dP1n_dt(kdx,1) = -((Delx_n^2)*as_n*Jni_des(kdx)*(x_sc^2)/(phisn_sc*sigmas_n*(epss_n^(Brugg_n)))) + ((P1_n(kdx+1) - 2*P1_n(kdx) + P1_n(kdx-1)));
    end
    
    dP1n_dt(N_n) = (P1_n(N_n) - P1_n(N_n -1)) + ((Delx_n)*I_app*x_sc/(S_n*phisn_sc*sigmas_n*(epss_n^(Brugg_n))))  ;
    
    %% Derivatives of electrolyte phase potential in positive electrode
    
    dP2p_dt(1) = P2_p(2) - P2_p(1);
    
    for kdx = 2:N_p-1
        dP2p_dt(kdx,1) = ((P2_p(kdx+1) - 2*P2_p(kdx) + P2_p(kdx-1))) + ((2*R*T*(transf - 1)*(log(C2_p(kdx+1)) - 2*log(C2_p(kdx)) +log(C2_p(kdx-1))))/(phie_sc*F)) ...
            + (((Delx_p^2)*as_p*Jpi(kdx)*(x_sc^2))/(phie_sc*kappae_p*(epse_p^(Brugg_p))));
    end
    
    dP2p_dt(N_p) =  ((kappae_m*(epse_m^(Brugg_m))*Delx_p)/(kappae_p*(epse_p^(Brugg_p))*Delx_m))*(P2_m(1) - P2_p(N_p)) - (P2_p(N_p) - P2_p(N_p - 1));
    
    %% Derivatives of electrolyte phase potential in membrane separator
    
    dP2m_dt(1) = (P2_m(2) - 2*P2_m(1) + P2_p(N_n)) + ((2*R*T*(transf - 1)*(log(C2_m(2)) - 2*log(C2_m(1)) +log(C2_p(N_p))))/(phie_sc*F)) ;
    
    for kdx = 2:N_m-1
        dP2m_dt(kdx,1) = (P2_m(kdx+1) - 2*P2_m(kdx) + P2_m(kdx-1)) + ((2*R*T*(transf - 1)*(log(C2_m(kdx+1)) - 2*log(C2_m(kdx)) +log(C2_m(kdx-1))))/(phie_sc*F));
    end
    
    dP2m_dt(N_m) = (P2_n(1) - 2*P2_m(N_m) + P2_m(N_m-1)) + ((2*R*T*(transf - 1)*(log(C2_n(1)) - 2*log(C2_m(N_m)) +log(C2_m(N_m-1))))/(phie_sc*F));
    
    %% Derivatives of electrolyte phase potential in negative electrode
    
    dP2n_dt(1) =  ((kappae_m*(epse_m^(Brugg_m))*Delx_n)/(kappae_n*(epse_n^(Brugg_n))*Delx_m))*(P2_n(1) - P2_m(N_m)) - (P2_n(2) - P2_n(1));
    
    if Index_chdisch == 1
        for kdx = 2:N_n-1
            dP2n_dt(kdx,1) = (P2_n(kdx+1) - 2*P2_n(kdx) + P2_n(kdx-1)) + ((2*R*T*(transf - 1)*(log(C2_n(kdx+1)) - 2*log(C2_n(kdx)) +log(C2_n(kdx-1))))/(phie_sc*F)) ...
                + (((Delx_n^2)*as_n*Jni_des(kdx)*(x_sc^2))/(phie_sc*kappae_n*(epse_n^(Brugg_n))));
        end
    elseif Index_chdisch == 4
        for kdx = 2:N_n-1
            dP2n_dt(kdx,1) = (P2_n(kdx+1) - 2*P2_n(kdx) + P2_n(kdx-1)) + ((2*R*T*(transf - 1)*(log(C2_n(kdx+1)) - 2*log(C2_n(kdx)) +log(C2_n(kdx-1))))/(phie_sc*F)) ...
                + (((Delx_n^2)*as_n*Jni_plating(kdx)*(x_sc^2))/(phie_sc*kappae_n*(epse_n^(Brugg_n))));
        end
    else
        for kdx = 2:N_n-1
            dP2n_dt(kdx,1) = (P2_n(kdx+1) - 2*P2_n(kdx) + P2_n(kdx-1)) + ((2*R*T*(transf - 1)*(log(C2_n(kdx+1)) - 2*log(C2_n(kdx)) +log(C2_n(kdx-1))))/(phie_sc*F)) ...
                + (((Delx_n^2)*as_n*Jni_sei(kdx)*(x_sc^2))/(phie_sc*kappae_n*(epse_n^(Brugg_n))));
        end
    end
    dP2n_dt(N_n) = P2_n(N_n);
    
    %% Derivatives of concentrations at solid electrolyte interphase
    
    dC3p_dt = C3_p - C1_p*Csp_sc/Csep_sc +(Jptot*Rs_p/(5*F*Ds_p*Csep_sc));
    
    if Index_chdisch == 1
        dC3n_dt = C3_n - C1_n*Csn_sc/Csen_sc +(Jntot*Rs_n/(5*F*Ds_n*Csen_sc));
    elseif Index_chdisch == 4
        dC3n_dt = C3_n - C1_n*Csn_sc/Csen_sc +(Jntot*Rs_n.*Jni0_des./(5*F*Ds_n*Csen_sc*Jni0_plating));
    else
        dC3n_dt = C3_n - C1_n*Csn_sc/Csen_sc +(Jntot*Rs_n.*Jni0_des./(5*F*Ds_n*Csen_sc*Jni0_sei));
    end
        
    %% Current density equations
    
    dJpidt = (Jpi - Jpi0.*(exp((alphaA*F*etap)/(R*T)) - exp(-((1-alphaA)*F*etap)/(R*T))))/abs(I_app/S_p);

    dJnidt_des = (Jni_des - Jni0_des.*(exp((alphaA*F*etan)/(R*T)) - exp(-((1-alphaA)*F*etan)/(R*T))))/abs(I_app/S_p);

    if Index_chdisch ~= 1
        dJnidt_sei = (Jni_sei - Jni0_sei.*(exp((alphaA*F*etan)/(R*T)) - exp(-((1-alphaA)*F*etan)/(R*T))))/abs(I_app/S_p);
    end
    if Index_chdisch == 4
        dJnidt_plating = (Jni_plating - Jni0_plating.*(exp((alphaA*F*etan)/(R*T)) - exp(-((1-alphaA)*F*etan)/(R*T))))/abs(I_app/S_p);
    end
    
    if Index_chdisch == 1
        
        %% Final expression
        dydt = [dC1p_dt; dC1n_dt; dC2p_dt; dC2m_dt; dC2n_dt; ...
            dP1p_dt; dP1n_dt; dP2p_dt; dP2m_dt; dP2n_dt;dC3p_dt;dC3n_dt;dJpidt;dJnidt_des];
        
    elseif Index_chdisch == 3
        
        %% voltage of cell
        dIapp_dt = P1_p(1)* phisp_sc - P1_n(end)*phisn_sc - Vcell_cc;
        
        %% Final expression
        dydt = [dC1p_dt; dC1n_dt; dC2p_dt; dC2m_dt; dC2n_dt; ...
            dP1p_dt; dP1n_dt; dP2p_dt; dP2m_dt; dP2n_dt; dC3p_dt; dC3n_dt; dJpidt; dJnidt_des; dJnidt_sei; dIapp_dt];
    elseif Index_chdisch == 4
        
        %% Final expression
        dydt = [dC1p_dt; dC1n_dt; dC2p_dt; dC2m_dt; dC2n_dt; ...
            dP1p_dt; dP1n_dt; dP2p_dt; dP2m_dt; dP2n_dt; dC3p_dt; dC3n_dt; dJpidt; dJnidt_des; dJnidt_sei; dJnidt_plating];

    else
        %% Final expression
        dydt = [dC1p_dt; dC1n_dt; dC2p_dt; dC2m_dt; dC2n_dt; ...
            dP1p_dt; dP1n_dt; dP2p_dt; dP2m_dt; dP2n_dt; dC3p_dt; dC3n_dt; dJpidt; dJnidt_des; dJnidt_sei];
        
    end
        
end

