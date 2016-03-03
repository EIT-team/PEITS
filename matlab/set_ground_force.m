function [gndpos] = set_ground_force(vtx,srf,grd_pos,PLOT)
% Usage: [gndpos] = set_ground_force(vtx,srf,grd_pos,cnts)
%
% General:
% Pick ground index surface between the 2 central electrodes
%
% Input: 
% vtx - vertices matrix {n x 3}
% srf - external surfaces {no. of external surfaces x 3}
% grd_pos - position of the ground node {1 x 3}
% PLOT - flag for plotting - 1 is on, 0 is off
%
% Output:
% gndpos - position of the ground vertex {1 x 3}
%-------------------------------------------------------------------------------------------------------------

srf_vtx = vtx(unique(reshape(srf,size(srf,1)*size(srf,2),1)),:);

grd_ind = grd_pos;
grd_ind_v = grd_ind (ones(1,size(srf_vtx,1)),:); % duplicate the electrode position, using tony's trick
dist_el = sum(((srf_vtx-grd_ind_v).^2)');
[~, grd_ind] = min(dist_el);
gndpos = srf_vtx(grd_ind,:);

if PLOT
%     figure, title('Mesh Surface');
%     trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',0.5);
%     colormap([0 0 0]);
%     daspect([1 1 1]);
%     
%     hold on;
%     axis image;
%     set(gcf,'Colormap',[0.6 0.7 0.9]);
%     
%     grid off
%     view(3);
    hold on;
    
    % plot on top of the surface mesh
    plot3(gndpos(1),gndpos(2),gndpos(3),'.', ...
        'MarkerSize',32,'Color','b');

    text(gndpos(1),gndpos(2),gndpos(3),'G','FontWeight','Bold','FontSize',16,'Color',[0.7 0.2 0.25]); % plots numbers on electrode pos
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This might become a part of the EIDORS suite
% Copyright (c) Lior Horesh 2004, EIT group, Medical Physics and Bioengineering, UCL, UK
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details
% EIDORS 3D version 2
% MATLAB Version 6.1.0.450 (R12.1) on PCWIN
% MATLAB License Number: 111672
% Operating System: Microsoft Windows Server 2003 Standard Edition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%