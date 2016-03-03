function [srf] = dubs3_2(tri)
% Usage: [srf] = dubs3_2(tri);
%
% General:
% Auxiliary function that extract the boundary faces of a given 3D volume (faster than the original dubs3.m)
%
% Input: 
% tri - simplices matrix {k x 4}
%
% Output:
% srf - outer boundary surfaces {no. of external surfaces x 3}
%------------------------------------------------------------------------------------------------------------------------

% create a list of all element facetts & sort
S = [tri(:,[1 2 3]); tri(:,[1 2 4]); tri(:,[2 3 4]); tri(:,[1 3 4])];
clear tri
N = sort(S,2);
clear S;
M = sortrows(N);
clear N;

i = 1;
inc = 1;

while i < size(M,1);
    ithrow = M(i,:);
    jthrow = M(i+1,:); 
    if ithrow == jthrow 
        i = i + 2;
    else
        % put in boundary node list
        srf(inc,:) = M(i,:);
        inc = inc + 1;
        i = i + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This might become a part of the EIDORS suite
% Copyright (c) Rebecca Yerworth & Lior Horesh 2004, EIT group, Medical Physics and Bioengineering, UCL, UK
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details
% EIDORS 3D version 2
% MATLAB Version 6.1.0.450 (R12.1) on PCWIN
% MATLAB License Number: 111672
% Operating System: Microsoft Windows Server 2003 Standard Edition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%