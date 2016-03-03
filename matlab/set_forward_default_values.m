function forward_settings = set_forward_default_values()

forward_settings.path = '/home/username/PEITS/dune-peits/';

% Mesh settings
forward_settings.mesh.name              = 'TA052_meters';
forward_settings.mesh.conductivities    = 0; % One value means uniform conductivity. 
                                             % Zero means values in the mesh. 
                                             % A vector encodes the different layers.
                                             % A string means that this conductivity file used
forward_settings.mesh.do_partitions          = false;

% Experimental settings
forward_settings.protocol           = 'current_protocol.txt'; % either the name of the file or a (n x 4) matrix
forward_settings.current            = 133e-6;   % current level in Ampere
forward_settings.contact_impedance  = 1e3;      % contact impedance of the electrodes
forward_settings.electrode_diameter = 7.0;      % diameter of the electrodes in milimeters

% Output
forward_settings.do_elec_volts      = true;     % write voltages to binary file
forward_settings.measORall          = true;     % write only the measured voltages or all electrode potentials
forward_settings.do_jacobian        = true;     % write jacobian matrix

% Perturbation
forward_settings.perturbation.do            = false;    % true for simulating a perturbation
forward_settings.perturbation.multORabs     = true;     % true for multiplication of original conductivity, false for absolute value
forward_settings.perturbation.value         = 1.5;      % either multiplied to original conductivity or set as absolute value
forward_settings.perturbation.radius        = 10;       % radius of perturbation in milimeters
forward_settings.perturbation.center        = [0.0, 0.0, 0.0]; % center of perturbation in meters
