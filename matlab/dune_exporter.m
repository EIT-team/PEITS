function [] = dune_exporter(Nodes, Tetra, sigma, filepath, filename, electrodepositions, groundposition)

% FOR EXPORTING IT TO MESH STRUCTURE RUN:
% dune_exporter(Nodes,Tetra,sigma,filepath,filename);
%
% filepath = '/home/username/PEITS/dune-peits/data/';
% filename = 'bla.dgf';
% electrodepositions (n x 3) - has default value for Markus' Mesh
% groundposition (1 x 3) - has default value for Markus' Mesh
%
% This function also determines the ground position precisely and writes a
% paramFILENAME file, which is needed by the dune forward solver
% Also, an electrode_positions_FILENAME file will be written.
%
% All files need to be placed into the 'dune-peits/data/' folder
%
% Author: Markus Jehl
%-----------------------------------------
if ~exist('groundposition','var')
    groundposition = [133.7000, 12.5200 55.0700]/1000;
end

if ~exist('electrodepositions','var')
    load('electrodepositions.mat');
end

if (size(Nodes,2)>3)
    Nodes = Nodes(:,1:3);
end

% Clean up the Nodes and Tetrahedra
todelete = length(Tetra)-length(find(Tetra==0))/4;
Tetra = Tetra(1:todelete,:); %cut away zeros at the end of Tetra
sigma = sigma(1:todelete); %the same for the conductivities
[Nodes, Tetra] = removeisolatednode(Nodes, Tetra);

% Get the ground position
srf = dubs3_2(Tetra);
elecpos = set_electrodes(Nodes, srf, electrodepositions, 1);
gndpos = set_ground_force(Nodes, srf, groundposition, 1);

% Find if conductivity is already correct or needs interpretation by forward
% solver (if all are integers)
if isinteger(sigma)
    isInt = 1;
else
    isInt = 0;
    
    % Check for case where sigma values are all whole numbers, but stored
    % as a double rather than an int.
    if (isempty(find(mod(sigma,1),1)))
        disp('INFO: All conductivity values are whole numbers.')
        disp('If these are supposed to be material reference (mat_refs) rather than absolute conducitivty values, convet sigma to integer representation')
        disp('(e.g. int8(sigma)) before passing to function')
    end
end

format longg
% Create file
inpfile = fopen([filepath filename],'w');

% Write Heading and Nodes
fprintf(inpfile,'DGF\n');
fprintf(inpfile,'vertex\n');
fprintf(inpfile,'firstindex 1\n');
fprintf(inpfile,'%6.18f %6.18f %6.18f  # %d\n', [Nodes' ; 1:size(Nodes,1)]); % write nodes and indices
fprintf(inpfile,'#\n');

% Write Tetrahedra
fprintf(inpfile,'Simplex\n');
fprintf(inpfile,'parameters 2\n');
fprintf(inpfile,'%d %d %d %d %f %d\n', [Tetra' ; sigma' ; 1:size(Tetra,1)]); % write tetrahedra and indices
fprintf(inpfile,'#\n');

fclose(inpfile);

% Write Electrode Position File
elecposfile = fopen([filepath 'electrode_positions_' filename(1:end-4)],'w');
fprintf(elecposfile,'%6.18f,%6.18f,%6.18f\n',elecpos');
fclose(elecposfile);

% Write Parameter File
paramfile = fopen([filepath 'param_' filename(1:end-4)],'w');

fprintf(paramfile, ['##### Parameters for ' filename(1:end-4) ' #####\n']);
fprintf(paramfile, ['# no. elements: ' num2str(size(Tetra,1)) ' # no. nodes: ' num2str(size(Nodes,1)) ' #\n\n']);

fprintf(paramfile, ['fem.io.macroGrid: ' filename '\n\n']);

fprintf(paramfile, ['electrode.use_node_assignment: false\n\n']);

fprintf(paramfile, ['electrode.positions: ../data/electrode_positions_' filename(1:end-4) '\n\n']);  % You might have to adjust this line later if you use different electrode positions!!
fprintf(paramfile, ['electrode.nodes: ../data/electrode_nodes_TA052_meters']); % If electrodes should be defined by nodes, change this file
fprintf(paramfile, ['surface.coordinates: ../data/surface_coordinates_TA052_precise']);

fprintf(paramfile, 'ground.hsquared: 1.5e-5\n');  % Also here you can let go of your creativity.
fprintf(paramfile, ['groundposition.x: ' num2str(gndpos(1), '%6.18f') '\n']);
fprintf(paramfile, ['groundposition.y: ' num2str(gndpos(2), '%6.18f') '\n']);
fprintf(paramfile, ['groundposition.z: ' num2str(gndpos(3), '%6.18f') '\n\n']);

if (isInt)
    fprintf(paramfile, 'fem.assign_conductivities: true');
else
    fprintf(paramfile, 'fem.assign_conductivities: false');
end
fclose(paramfile);
