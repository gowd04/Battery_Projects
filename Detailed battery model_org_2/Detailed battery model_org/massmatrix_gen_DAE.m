function M  = massmatrix_gen_DAE(~,~,Index_chdisch,~,~,~,~,~,~,~)

    parameters2
    
    if Index_chdisch == 1
        M = diag(ones(6*(N_p+N_n) + 2*N_m ,1));
    else
        M = diag(ones(6*(N_p+N_n) + 2*N_m +N_n,1));
    end
    
    % Derivatives of electrolyte phase conc in positive electrode
    M(N_p + N_n+1, N_p + N_n+1)     = 0;
    M(2*N_p + N_n, 2*N_p + N_n)     = 0;
    
    % Derivatives of electrolyte phase conc in membrane separator
    % all diagonal elements will be 1. so no change
    
    % Derivatives of electrolyte phase conc in negative electrode
    M(2*N_p + N_m + N_n + 1, 2*N_p + N_m + N_n + 1)         = 0;
    M(2*N_p + N_m + 2*N_n , 2*N_p + N_m + 2*N_n )           = 0;
    
    % Derivatives of solid phase potential in positive and negative electrode    
    for jdx = (2*N_p + N_m + 2*N_n + 1) :(3*N_p + N_m + 3*N_n)
        M(jdx,jdx) = 0;
    end
    
    % Derivatives of electrolyte phase potential    
    for jdx = 3*N_p + N_m + 3*N_n + 1: 4*N_p + 2*N_m + 4*N_n
        M(jdx,jdx) = 0;
    end
    
    % Derivatives of concentration at solid electrolyte interphase
    for jdx = 4*N_p + 2*N_m + 4*N_n + 1: 5*N_p + 2*N_m + 5*N_n
        M(jdx,jdx) = 0;
    end
    
    % Current density equations    
    
    if Index_chdisch == 1
        for jdx = 5*N_p + 2*N_m + 5*N_n+1:6*N_p + 2*N_m + 6*N_n
            M(jdx,jdx) = 0;
        end
    elseif Index_chdisch == 4
        for jdx = 5*N_p + 2*N_m + 5*N_n+1:6*N_p + 2*N_m + 8*N_n
            M(jdx,jdx) = 0;
        end    
    else
        for jdx = 5*N_p + 2*N_m + 5*N_n+1:6*N_p + 2*N_m + 7*N_n
            M(jdx,jdx) = 0;
        end        
    end
  
    if Index_chdisch == 3
        M(6*N_p + 2*N_m + 7*N_n+1,6*N_p + 2*N_m + 7*N_n+1) = 0;
    end
end











