function [id, sigma] = load_sigma_vector_binary(filename)

  fid = fopen(filename,'r');
  % if standard reading is not the correct format for a given binary
  % file, activate the following:
  %fid = fopen(filename,'r','ieee-be');
  
  magicstr = char(fread(fid,3,'char'))';
  if ~isequal(magicstr,'DDV')
    error('read magicstr doe not indicate Dune Dof Vector!');
  end;
  
  magicint = fread(fid,1,'int');
  magicdouble = fread(fid,1,'double');
  
  if (magicint~=111) | (magicdouble~=111.0)
    error(['magic numbers not read correctly, change the binary format in' ...
	   ' this reading routine!']);
  end;
  
  nentries = fread(fid,1,'int');
  
  disp(['generating and reading vector with ',num2str(nentries),' entries.']);

  v = fread(fid,2*nentries,'double');

  eofstr = char(fread(fid,3,'char'))';
  if ~isequal(eofstr,'EOF')
    error('read eofstr does not indicate end of binary file!');
  end;
  
  id = v(1:2:end);
  sigma = v(2:2:end);

  fclose(fid);
