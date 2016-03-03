function [v,j,output] = run_forward_solver(forward_settings, processes)

% write the settings to the parameter files
if ~exist([forward_settings.path, 'data/', forward_settings.mesh.name, '.dgf'], 'file')
    disp('This mesh has not yet been written to DGF format. Use dune_exporter.m to do that.');
    return
end
if ~exist([forward_settings.path, 'partitions/', forward_settings.mesh.name, '.dgf.', num2str(processes), '_0'], 'file')
    disp('This mesh has not yet been partitioned to the defined number of processes. This will take a while.');
    forward_settings.mesh.do_partitions = true;
end
if isa(forward_settings.protocol, 'double')
    disp('Creating new protocol file "data/protocol.txt".');
    fid = fopen([forward_settings.path,'data/protocol.txt'],'w');
    fprintf(fid,'%f,%f,%f,%f\n',forward_settings.protocol');
    fclose(fid);
    forward_settings.protocol = 'protocol.txt';
end

write_parameter_files(forward_settings);

% run the solver
[status,output] = system(['mpirun -np ', num2str(processes), ' ', forward_settings.path, 'src/dune_peits']);

if status
    disp('Something went wrong when running the solver:');
    disp(output);
end

% read out the results
if forward_settings.do_elec_volts
    % find out which is the newest file and load it
    d = dir([forward_settings.path,'output/electrodevoltages*.bin']);
    [~,index] = max([d.datenum]);
    v = load_electrode_voltages_binary([forward_settings.path,'output/',d(index).name]);
else
    v = 0;
end

if forward_settings.do_jacobian
    % find out which is the newest file and load it. also, sort the
    % jacobian by the element ID's
    d = dir([forward_settings.path,'output/sigmavector*.bin']);
    [~,index] = max([d.datenum]);
    [id,~] = load_sigma_vector_binary([forward_settings.path,'output/',d(index).name]);
    d = dir([forward_settings.path,'output/jacobian*.bin']);
    [~,index] = max([d.datenum]);
    j_unsorted = load_jacobian_binary([forward_settings.path,'output/',d(index).name]);
    j = zeros(size(j_unsorted));
    j(:,id) = j_unsorted;
else
    j = 0;
end
