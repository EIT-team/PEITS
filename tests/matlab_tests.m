addpath ('../matlab');

% Set up
load('../resources/TA052.mat') % Load mesh
sigma_values = mat_ref;
filepath = './';
filename = 'test.out';


%% Tests
% Expected values
fem_true = 'fem.assign_conductivities: true'
fem_false = 'fem.assign_conductivities: false'

disp('Check that doubles are interpreted as conducitivity values and integers as mat_refs')

% Test 1
test = 'doubles'; % Pass in sigma as doubles
disp(['Testing: ' test])
dune_exporter(vtx, tri, sigma_values, filepath, filename, pos)
file = textread('param_test', '%s', 'delimiter', '\n');
num_lines = size(file,1);

line_to_test = file{num_lines};
assert( isequal(line_to_test,fem_false),['Error when testing: ' test])
assert( ~isequal(line_to_test,fem_true),['Error when testing: ' test])

% Test 2
test = 'whole numbers as doubles'; % Pass in sigma as integer values, but stored as a double
disp(['Testing: ' test])

sigma_values = ceil(sigma_values);
dune_exporter(vtx, tri, sigma_values, filepath, filename, pos)
file = textread('param_test', '%s', 'delimiter', '\n');
num_lines = size(file,1);

line_to_test = file{num_lines};
assert( isequal(line_to_test,fem_false),['Error when testing: ' test])
assert( ~isequal(line_to_test,fem_true),['Error when testing: ' test])

% Test 3
test = 'integers'; % Pass in integers
disp(['Testing: ' test])
sigma_values = int8(sigma_values);
dune_exporter(vtx, tri, sigma_values, filepath, filename, pos)
file = textread('param_test', '%s', 'delimiter', '\n');
num_lines = size(file,1);

line_to_test = file{num_lines};

assert( isequal(line_to_test,fem_true),['Error when testing: ' test])
assert( ~isequal(line_to_test,fem_false),['Error when testing: ' test])


% Clean up
fclose('all'); % Close open files
clear fem* file* line_to_test mat_ref mesh* num_lines pos sigma_values test tri vtx % Remove variables
delete('param_test', 'test.out', 'electrode_positions_test')