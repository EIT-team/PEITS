function [coordinates, origin] = generate_surface_coords (Mesh)
%% Mesh has to have faces defined. we use dubs3_2 for doing that
vtx = Mesh.Nodes_faces;
srf = Mesh.Faces;

[max_z,ind_max_z] = max(vtx(:,3));
min_z = min(vtx(:,3));

% we only use the upper 2/3 of the head (z-coord).
%vtx=vtx(vtx(:,3)>=min_z+(max_z-min_z)/3,:);
%[max_z,ind_max_z] = max(vtx(:,3));
vtx=vtx(vtx(:,3)>=min_z+(max_z-min_z)/6,:);
[max_z,ind_max_z] = max(vtx(:,3));

%%
% define the distance between coordinate system points
point_distance=8;
eps = 2; % precision of coordinate points

% define the center of the surface coordinate system as the one with largest z coordinate
coord0 = vtx(ind_max_z,:);

% Select first stripe parallel to x axis
sel=vtx(abs(vtx(:,2)-coord0(2))<eps,:);

% then we remove all vertices close to the origin since we do not consider them later
dist=(sum((sel-repmat(coord0,length(sel),1)).^2,2)).^0.5;
sel(dist<point_distance-eps,:)=[];

coordinates(1,:)=coord0;   % this coordinate matrix is going to be extended in all directions with new coordinate points
origin_m = 1; origin_n = 1; % m is the x position of the point in the coordinates matrix and n the y position
[el_f] = make_stripe(coord0, 1, 2, sel, eps, point_distance); % stripe of holes along x
[el_b] = make_stripe(coord0, -1, 2, sel, eps, point_distance); % stripe of holes along x
coordinates=[el_b;coordinates;el_f];
origin_m = origin_m + size(el_b,1);

% Select the first stripe parallel to y
sel=vtx(abs(vtx(:,1)-coord0(1))<eps,:);

% remove elements from queue
dist=(sum((sel-repmat(coord0,length(sel),1)).^2,2)).^0.5;
sel(dist<point_distance-eps,:)=[];

[el_f] = make_stripe(coord0, 2, 1, sel, eps, point_distance); % stripe of holes along y forward
[el_b] = make_stripe(coord0, -2, 1, sel, eps, point_distance); % stripe of holes along y backward
tmp_matrix = inf*ones(size(coordinates,1),3*size(el_f,1)+3*size(el_b,1));
tmp_matrix(origin_m,:) = [reshape(el_b',1,[]),reshape(el_f',1,[])];
coordinates = [tmp_matrix(:,1:3*size(el_b,1)),coordinates,tmp_matrix(:,3*size(el_b,1)+1:end)];
origin_n = origin_n + 3*size(el_b,1);

 
%% now filling the quarters
sel = vtx;
coordinates = fill_quarters(coordinates, origin_m, origin_n, sel, eps, point_distance);


%% Plotting
p=1:size(coordinates,1);
for j=1:size(coordinates,2)/3
    scatter3(0.001*coordinates(p,(j-1)*3+1),0.001*coordinates(p,(j-1)*3+2),0.001*coordinates(p,(j-1)*3+3),'r');
    hold on;
end
daspect([1 1 1]);

%% output
origin = [origin_m, origin_n];

end

%%
function [points] = make_stripe(first, dim_along, dim_cross, sel, eps, point_distance)
% dim_along: along which dimension; sign: negative=backwards, positive=forwards;
% first: first coordinate (not included in output)
% sel: the stripe of candidate vertices
points(1,:)=first;
dd=abs(dim_along);

next=first;
status=0;
while (status~=2)
    dist=(sum((sel-repmat(next,size(sel,1),1)).^2,2)).^0.5;
    
    if (status==0)     % place the first electrode in the direction
        if dim_along>0 % direction of placement
            p=find(abs(dist-point_distance)<eps & sel(:,dd)>next(dd)); % candidates
            p1=find(dist<2*point_distance-eps | sel(:,dd)<next(dd));   % delete these
        else
            p=find(abs(dist-point_distance)<eps & sel(:,dd)<next(dd));
            p1=find(dist<2*point_distance-eps | sel(:,dd)>next(dd));
        end
        status=1;
    else % all other electrodes
        p=find(abs(dist-point_distance)<eps);
        p1=find(dist<2*point_distance-eps);
    end
    if (~isnan(p)) % continue as long as there are candidates
        next=sel(p,:);
        sel(p1,:)=[];
        next=mean(next,1);
        next(dim_cross)=first(dim_cross);
        points = [points;next];
    else
        status=2;       
    end
end
points(1,:)=[]; % we do not need first here
if dim_along<0
    points = flipud(points);
end

end

%% Function to fill all four quarters
function coordinates = fill_quarters (coordinates, origin_m, origin_n, sel, eps, point_distance)
% this function fills the four quarters spanned by the x and y coordinates
% given in 'coordinates'
x_axis = reshape(coordinates(origin_m,:),3,[])';
y_axis = coordinates(:,origin_n:origin_n+2);
origin = coordinates(origin_m,origin_n:origin_n+2);

all_coords = [x_axis;y_axis];
for i=1:size(all_coords,1) % remove all non-candidates from selected vertices
    dist=(sum((sel-repmat(all_coords(i,:),length(sel),1)).^2,2)).^0.5;
    p=find(dist<point_distance-eps);
    sel(p,:)=[];
end

% separate the vertex selection into the four quarters
sel_p = sel(sel(:,1)>origin(1),:);
sel_pp = sel_p(sel_p(:,2)>origin(2),:);
sel_pm = sel_p(sel_p(:,2)<origin(2),:);
sel_m = sel(sel(:,1)<origin(1),:);
sel_mp = sel_m(sel_m(:,2)>origin(2),:);
sel_mm = sel_m(sel_m(:,2)<origin(2),:);

% fill x+ y+ quarter
coordinates = fill_quarter(1,1,x_axis((origin_n+2)/3:end,:),y_axis(origin_m:end,:),coordinates,origin_m,origin_n,sel_pp,eps,point_distance);

% fill x+ y- quarter
coordinates = fill_quarter(1,-1,x_axis((origin_n+2)/3:end,:),y_axis(1:origin_m,:),coordinates,origin_m,origin_n,sel_mp,eps,point_distance);

% fill x- y+ quarter
coordinates = fill_quarter(-1,1,x_axis(1:(origin_n+2)/3,:),y_axis(origin_m:end,:),coordinates,origin_m,origin_n,sel_pm,eps,point_distance);

% fill x- y- quarter
coordinates = fill_quarter(-1,-1,x_axis(1:(origin_n+2)/3,:),y_axis(1:origin_m,:),coordinates,origin_m,origin_n,sel_mm,eps,point_distance);

end

%% Function to fill a specific quarter
function coordinates = fill_quarter(x_dir, y_dir, x_axis, y_axis, coordinates, origin_m, origin_n, sel, eps, point_distance)

if x_dir>0
    x_along = x_axis(2:end,:);
else
    x_along = flipud(x_axis(1:end-1,:));
end

if y_dir>0
    y_along = y_axis(2:end,:);
else
    y_along = flipud(y_axis(1:end-1,:));
end

for i=1:size(x_along,1)
    newline_x = inf*ones(size(y_along,1)+1,3);
    newline_x(1,:) = x_along(i,:);
    for j=1:size(y_along,1)
        
        dist_x=(sum((sel-repmat(newline_x(j,:),size(sel,1),1)).^2,2)).^0.5;
        dist_y=(sum((sel-repmat(y_along(j,:),size(sel,1),1)).^2,2)).^0.5;
        
        p=find(abs(dist_x-point_distance)<eps & abs(dist_y-point_distance)<eps); % find the intersection of the two stripes
       
        if (~isnan(p))
            next=sel(p,:);
            next=mean(next,1);
            
            newline_x(j+1,:)=next;
            dist=(sum((sel-repmat(next,size(sel,1),1)).^2,2)).^0.5;
            sel(dist<=point_distance-eps,:)=[];
        else
            break
        end
        
    end
    
    y_along = newline_x(2:end,:);
    
    % write the new line into the coordinate matrix
    if y_dir>0
        coordinates(origin_m+1:end,origin_n + x_dir*3*i:origin_n + x_dir*3*i + 2) = newline_x(2:end,:);
    else
        coordinates(1:origin_m-1,origin_n + x_dir*3*i:origin_n + x_dir*3*i + 2) = flipud(newline_x(2:end,:));
    end
end

end
