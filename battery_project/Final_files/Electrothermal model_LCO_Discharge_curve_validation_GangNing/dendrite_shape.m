% clc
% clear
% 
% load('dendrite.mat')
% indices = find(Li_den_matrix(:,4) ~= 0 & Li_den_matrix(:,1)>=0 & Li_den_matrix(:,1)<= 2.5e-5 & Li_den_matrix(:,2)>=0 & Li_den_matrix(:,2)<=2.5e-5 );

indices = find(Li_den_matrix(:,4) ~= 0 );%
max_len = max(Li_den_matrix(indices,4));
min_len = min(Li_den_matrix(indices,4));

% figure
% [x1,x2] = meshgrid(0:1e-4:0.001,0:1e-4:0.001);
% h = surf(x1,x2,zeros(size(x1)));
% set(h, 'FaceColor',[102 133 134]./ 255, 'FaceAlpha',0.5, 'EdgeAlpha', 0);
% for ndx = 1:length(indices)
%     jdx = indices(ndx);
%     Len = Li_den_matrix(jdx,4);
%     r   = Li_den_matrix(jdx,3);
%     
%     sc =  (Len - min_len)/(max_len - min_len)+eps;
%     
%     if Len > r
%         x0  = [Li_den_matrix(jdx,1:2) 0];
%         
%         n = 100; m = 200;
%         
%         theta1  = linspace(0, 2*pi, n)';
%         x1 = repmat(r*cos(theta1) + x0(1),1,m);
%         y1 = repmat(r*sin(theta1) + x0(2),1,m);
%         z1 = repmat(linspace(0,Len-r,m),n,1);
%         
%         theta1 = repmat(linspace(0, 2*pi, n)',1,m);
%         theta2 = repmat(linspace(0, pi/2, m),n,1);
%         x2 = r*cos(theta1).*sin(theta2) + x0(1);
%         y2 = r*sin(theta1).*sin(theta2) + x0(2);
%         z2 = r*cos(theta2) + x0(3) + (Len-r);
%         
%         hold all
%         h(1) = surf(x1,y1,z1);
%         
%         h(2) = surf(x2,y2,z2);
%         set(h, 'FaceColor',[sc*102 133 sc*134]./ 255, 'FaceAlpha',1, 'EdgeAlpha', 0);
%         
%     else
%         x0  = [Li_den_matrix(jdx,1:2) 0];
%         
%         n = 100; m = 2000;
%         
%         phi_1 = asin(r/sqrt(2*r^2 + Len^2 - 2*r*Len));
%         phi_2 = -asin(r/sqrt(2*r^2 + Len^2 - 2*r*Len));
%         
%         theta1 = repmat(linspace(0, 2*pi, n)',1,m);
%         theta2 = repmat(linspace(0 , pi/2, m),n,1);
%         x2 = r*cos(theta1).*sin(theta2) + x0(1);
%         y2 = r*sin(theta1).*sin(theta2) + x0(2);
%         z2 = r*cos(theta2) + x0(3) - (r-Len);
%         
%         ind = find(z2(1,:) >= 0);
%         
%         x2 = x2(:,ind);
%         y2 = y2(:,ind);
%         z2 = z2(:,ind);
%         
%         hold on
%         
%         h = surf(x2,y2,z2);
% 
%         set(h,'FaceColor',[sc*102 133 sc*134] ./ 255,'FaceAlpha',1, 'EdgeAlpha', 0);
%     end
% end
