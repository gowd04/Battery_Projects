function [value,isterminal,direction] = breakode_battery_new(t,y,Index_chdisch,~,~,~,~,Cs_max_p,Cs_max_n,~,~, EOCV, ~)

    parameters2_NMC
    
    C1_p_act = y(1:N_p)*Csp_sc/Cs_max_p;
    C1_n_act = y(N_p+1:N_p+N_n)*Csn_sc/Cs_max_n;
    
    P1_p = y(2*N_p+N_m+2*N_n+1:3*N_p+N_m+2*N_n);
    P1_n = y(3*N_p+N_m+2*N_n+1:3*N_p+N_m+3*N_n);
    
    
    if Index_chdisch == 1 %% CC discharging
        
        if (sign(P1_p(1) - P1_n(end) - EODV) + 1) == 0 || ~isempty(find(C1_p_act >= 0.9999, 1)) || ~isempty(find(C1_n_act <= 0.0001, 1))
            value = 0;
        else
            value = 1;
        end
        isterminal = 1;
        direction  = 0;
        
    elseif Index_chdisch == 4 %% CC charging
        
        if (sign(P1_p(1) - P1_n(end) - EOCV) - 1) == 0 || ~isempty(find(C1_n_act >= 0.9, 1)) || ~isempty(find(C1_p_act <= 0.394, 1))
            value = 0;
        else
            value = 1;
        end
        isterminal = 1;
        direction  = 0;
        
    elseif Index_chdisch == 3 %% CV charging
        
        I_app = y(7*N_p+3*N_m+8*N_n+1);
        if (sign(I_app - EOCC) + 1) == 0 || ~isempty(find(C1_n_act >= 0.9, 1)) || ~isempty(find(C1_p_act <= 0.394, 1))
            value = 0;
        else
            value = 1;
        end
        isterminal = 1;
        direction  = 0;
        
%     elseif Index_chdisch == 4 %% CC charging with Li plating
%         
%         if (sign(P1_p(1) - P1_n(end) - (EOCV + 0.2)) - 1) == 0 || ~isempty(find(C1_p_act >= 0.9999, 1)) || ~isempty(find(C1_n_act <= 0.0001, 1))
%             value = 0;
%         else
%             value = 1;
%         end
%         isterminal = 1;
%         direction  = 0;
        
    end
end

