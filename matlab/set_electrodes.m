function [elecpos] = set_electrodes(vtx, srf, pos, fig)
% Usage: [elec, cnts, srf_area] = set_electrodes(vtx, srf, pos, w, method, fig);
%
% General:
% Project electrodes on the mesh surface, starting from the element centre which is closest to the electrode position
%
% Input:
% vtx - vertices matrix {n x 3}
% srf - external surfaces {no. of external surfaces x 3}
% pos - position of the electrodes surfaces vertices {no. of electrodes x 3}
% fig - flag, if not exist surface mesh will be ploted {} (optional)
%
% Output:
% elecpos - centres of surface triangles closest to given electrode positions
%-------------------------------------------------------------------------------------------------------------

% vector containing the geometric center of each triangular  surface in x,y,z coordinates.
cnts = (vtx(srf(:,1),:) + vtx(srf(:,2),:) + vtx(srf(:,3),:))/3;
elecpos = zeros(size(pos));
elec_srf = zeros(size(pos,1),1);
for j = 1 : size(pos,1)
        
    % finds the closest surface to the electrode position
    el_pos = pos(j,:);
    el_pos_v = el_pos (ones(1,length(cnts)),:); % duplicate the electrode position, using tony's trick
    dist_el = sum(((cnts-el_pos_v).^2)');
    elec_srf(j) = find(dist_el == min(dist_el)); % finds the index of the closest surface to the electrode

    % write the centre of the closest surface triangle
    elecpos(j,:) = cnts(elec_srf(j),:);
end

%-------------------------------plot the electrodes---------------------------------------
if fig
    figure, title('electrodes');
    trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',0.5);
    colormap([0 0 0]);
    daspect([1 1 1]);
    
    hold on;
    axis image;
    set(gcf,'Colormap',[0.6 0.7 0.9]);
    
    grid off
    view(3);
    hold on
    
    l = srf(elec_srf,1); m = srf(elec_srf,2); n = srf(elec_srf,3);
    
    Xs = [vtx(l,1)';vtx(m,1)';vtx(n,1)'];
    Ys = [vtx(l,2)';vtx(m,2)';vtx(n,2)'];
    Zs = [vtx(l,3)';vtx(m,3)';vtx(n,3)'];
    
    patch(Xs,Ys,Zs,'b');
    
    % add numbers to the electrodes
    el_caption = 1:size(pos,1);
    hold on;
    for i = 1:size(pos,1)
        text(pos(i,1),pos(i,2),pos(i,3),num2str(el_caption(i)),'FontWeight','Bold','FontSize',16,'Color',[0.4 0.2 0.65]); % plots numbers on electrode pos
    end
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