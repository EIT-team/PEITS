function write_parameter_files(forward_settings)

% rename the original parameter file
movefile([forward_settings.path, 'data',filesep,'parameter'],[forward_settings.path, 'data',filesep,'parameter~']);

% open parameter files
oldfile = fopen([forward_settings.path, 'data',filesep,'parameter~'], 'r');
newfile = fopen([forward_settings.path, 'data',filesep,'parameter'], 'w');

% write changed settings into new parameter file
tline=fgetl(oldfile);
while ischar(tline)
    
    if strfind(tline, 'fem.io.loadPartitions:')
        if forward_settings.mesh.do_partitions
            fprintf(newfile, '%s\n', 'fem.io.loadPartitions: false');
        else
            fprintf(newfile, '%s\n', 'fem.io.loadPartitions: true');
        end
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'fem.io.do_elec_volts:')
        if forward_settings.do_elec_volts
            fprintf(newfile, '%s\n', 'fem.io.do_elec_volts: true');
        else
            fprintf(newfile, '%s\n', 'fem.io.do_elec_volts: false');
        end
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'fem.io.do_jacobian:')
        if forward_settings.do_jacobian
            fprintf(newfile, '%s\n', 'fem.io.do_jacobian: true');
        else
            fprintf(newfile, '%s\n', 'fem.io.do_jacobian: false');
        end
        tline=fgetl(oldfile);
        continue
    end
    
    if any(strfind(tline, 'mesh:') == 1)
        fprintf(newfile, '%s\n', ['mesh: ',forward_settings.mesh.name]);
        tline=fgetl(oldfile);
        continue
    end

    if strfind(tline, 'fem.io.load_sigma_file_separately:')
      if forward_settings.fem.io.load_sigma_separately
        fprintf(newfile, '%s\n', 'fem.io.load_sigma_separately: true');
      else
        fprintf(newfile, '%s\n', 'fem.io.load_sigma_separately: false');
      end
      tline=fgetl(oldfile);
      continue;
    end

    if strfind(tline, 'fem.io.separate_sigma_file:')
      if ~isempty(forward_settings.fem.io.separate_sigma_file)
        fprintf(newfile, '%s\n', ['fem.io.separate_sigma_file: ', forward_settings.fem.io.separate_sigma_file]);
      else
        fprintf(newfile, '%s\n', 'fem.io.separate_sigma_file: none');
      end
      tline=fgetl(oldfile);
      continue;
    end
    
    if strfind(tline, 'conductivities:')
        if (length(forward_settings.mesh.conductivities)>1)&&isa(forward_settings.mesh.conductivities, 'char')
            fprintf(newfile, '%s\n', ['conductivities: ..',filesep,'data',filesep,...
                    forward_settings.mesh.conductivities]);
            tline=fgetl(oldfile);
            continue
        else if (length(forward_settings.mesh.conductivities)>1)&&isa(forward_settings.mesh.conductivities, 'double')
            disp('Writing new conductivity file "data/conductivities".');
            fid = fopen([forward_settings.path,'data',filesep,'conductivities'],'w');
            fprintf(fid,'%f\n',forward_settings.mesh.conductivities);
            fclose(fid);
            fprintf(newfile, '%s\n', 'conductivities: ..',filesep,'data',filesep,'conductivities');
            tline=fgetl(oldfile);
            continue
            end
        end
    end
    
    if any(strfind(tline, 'current.protocol:') == 1)
        fprintf(newfile, '%s\n', ['current.protocol: ..',filesep,'data',filesep,forward_settings.protocol]);
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'mesh.perturbation:')
        if forward_settings.perturbation.do
            fprintf(newfile, '%s\n', 'mesh.perturbation: true');
        else
            fprintf(newfile, '%s\n', 'mesh.perturbation: false');
        end
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'mesh.perturbation.multORabs:')
        if forward_settings.perturbation.multORabs
            fprintf(newfile, '%s\n', 'mesh.perturbation.multORabs: true');
        else
            fprintf(newfile, '%s\n', 'mesh.perturbation.multORabs: false');
        end
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'mesh.perturbation.value:')
        fprintf(newfile, '%s\n', ['mesh.perturbation.value: ',num2str(forward_settings.perturbation.value)]);
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'mesh.perturbation.radius:')
        fprintf(newfile, '%s\n', ['mesh.perturbation.radius: ',num2str(forward_settings.perturbation.radius)]);
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'mesh.perturbation.pos_x:')
        fprintf(newfile, '%s\n', ['mesh.perturbation.pos_x: ',num2str(forward_settings.perturbation.center(1))]);
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'mesh.perturbation.pos_y:')
        fprintf(newfile, '%s\n', ['mesh.perturbation.pos_y: ',num2str(forward_settings.perturbation.center(2))]);
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'mesh.perturbation.pos_z:')
        fprintf(newfile, '%s\n', ['mesh.perturbation.pos_z: ',num2str(forward_settings.perturbation.center(3))]);
        tline=fgetl(oldfile);
        continue
    end
    
	fprintf(newfile, '%s\n', tline);
    tline=fgetl(oldfile);
end

% close both files
fclose(oldfile);
fclose(newfile);


% rename the original standardparams file
movefile([forward_settings.path, 'data',filesep,'standardparams'],[forward_settings.path, 'data',filesep,'standardparams~']);

% open parameter files
oldfile = fopen([forward_settings.path, 'data',filesep,'standardparams~'], 'r');
newfile = fopen([forward_settings.path, 'data',filesep,'standardparams'], 'w');

% write changed settings into new parameter file
tline=fgetl(oldfile);
while ischar(tline)
    
    if strfind(tline, 'fem.io.write_only_measured_voltage:')
        if forward_settings.measORall
            fprintf(newfile, '%s\n', 'fem.io.write_only_measured_voltage: true');
        else
            fprintf(newfile, '%s\n', 'fem.io.write_only_measured_voltage: false');
        end
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'fem.uniform_conductivity:')
        if (length(forward_settings.mesh.conductivities)==1)&&forward_settings.mesh.conductivities
            fprintf(newfile, '%s\n', 'fem.uniform_conductivity: true');
        else
            fprintf(newfile, '%s\n', 'fem.uniform_conductivity: false');
        end
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'fem.uniform_conductivity_value:')
        if (length(forward_settings.mesh.conductivities)==1)
            fprintf(newfile, '%s\n', ['fem.uniform_conductivity_value: ',num2str(forward_settings.mesh.conductivities)]);
            tline=fgetl(oldfile);
            continue
        end
    end
    
    if strfind(tline, 'contact.impedance:')
        fprintf(newfile, '%s\n', ['contact.impedance: ',num2str(forward_settings.contact_impedance)]);
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'input.current:')
        fprintf(newfile, '%s\n', ['input.current: ',num2str(forward_settings.current)]);
        tline=fgetl(oldfile);
        continue
    end
    
    if strfind(tline, 'electrode.diameter:')
        fprintf(newfile, '%s\n', ['electrode.diameter: ',num2str(forward_settings.electrode_diameter)]);
        tline=fgetl(oldfile);
        continue
    end
    
	fprintf(newfile, '%s\n', tline);
    tline=fgetl(oldfile);
end

% close both files
fclose(oldfile);
fclose(newfile);
