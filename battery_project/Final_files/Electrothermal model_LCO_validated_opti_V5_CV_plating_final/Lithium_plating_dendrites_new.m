function Li_den_matrix = Lithium_plating_dendrites_new(Li_den_matrix, Li_plated_volume)
    
    V_rem   = Li_plated_volume;
    N_sites = size(Li_den_matrix,1);
    len_max1 = 5e-6; % Maximum possible length of dendrite during one deposition
    
    while V_rem > 0
        ind_site = randi(N_sites);
        l_den = rand(1)*len_max1;
        r_den = Li_den_matrix(ind_site,3);
        if Li_den_matrix(ind_site,4) + l_den < r_den
            V_den = (pi*(Li_den_matrix(ind_site,4) + l_den)*(3*r_den^2  +  (Li_den_matrix(ind_site,4) + l_den)^2)/6)  - Li_den_matrix(ind_site,5);
            V_rem_New = V_rem - V_den;
            
            if V_rem_New < 0
                l_den = fsolve(@(x) V_rem - ((pi*(Li_den_matrix(ind_site,4) + x)*(3*r_den^2  +  (Li_den_matrix(ind_site,4) + x)^2)/6)  - Li_den_matrix(ind_site,5)),0,optimoptions('fsolve','Display','None','TolFun',1e-25));
                V_den = (pi*(Li_den_matrix(ind_site,4) + l_den)*(3*r_den^2  +  (Li_den_matrix(ind_site,4) + l_den)^2)/6)  - Li_den_matrix(ind_site,5);
                
                Li_den_matrix(ind_site,4) = Li_den_matrix(ind_site,4) + l_den;
                Li_den_matrix(ind_site,5) = Li_den_matrix(ind_site,5) + V_den;
                
                break
            end
            
        else
            V_den = (2/3)*pi*r_den^3  +  pi*(r_den^2)*(l_den - r_den);
            V_rem_New = V_rem - V_den;
            
            if V_rem_New < 0
                l_den = ((V_rem - (2/3)*pi*r_den^3)/( pi*(r_den^2))) + r_den;
                V_den = (2/3)*pi*r_den^3  +  pi*(r_den^2)*(l_den - r_den);
                
                Li_den_matrix(ind_site,4) = Li_den_matrix(ind_site,4) + l_den;
                Li_den_matrix(ind_site,5) = Li_den_matrix(ind_site,5) + V_den;
                
                break
            end
        end        
        
        V_rem = V_rem - V_den;
        
        Li_den_matrix(ind_site,4) = Li_den_matrix(ind_site,4) + l_den;
        Li_den_matrix(ind_site,5) = Li_den_matrix(ind_site,5) + V_den;
        
    end
end

%% Generation of possible sites
% B_height = 0.878; % Battery height
% B_depth  = 0.1; % Battery depth
% r_max    = 1e-4; % Maximum possible radius of Li dendrite
% N_max    = 5000; % Maximum number of possible sites on electrode surface
% r_min    = r_max/20; % Maximum possible radius of Li dendrite

% [x_all, y_all, r_all] = possible_sites(B_depth, B_height, r_max,r_min, N_max+1);
% Li_den_matrix = [x_all, y_all, r_all];
% Li_den_matrix(:,4:5) = 0;
% clearvars -except Li_den_matrix

% close all

% for jdx = 1:length(Li_den_matrix(:,1))
%     circle([Li_den_matrix(jdx,1), Li_den_matrix(jdx,2)],Li_den_matrix(jdx,3))
%     hold all
% end
% 

% B_height = 0.001; % Battery height
% B_depth  = 0.001; % Battery depth
% r_max    = 1e-4; % Maximum possible radius of Li dendrite
% r_min    = r_max/10; % Maximum possible radius of Li dendrite
% 
% N_max    = 20; % Maximum number of possible sites on electrode surface
% [x_all, y_all, r_all] = possible_sites(B_depth, B_height, r_max, r_min,N_max+1);
% Li_den_matrix = [x_all, y_all, r_all];
% Li_den_matrix(:,4:5) = 0;
% clearvars -except Li_den_matrix

