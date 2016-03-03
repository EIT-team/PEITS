function A = load_electrode_voltages_binary(filename)

  fid = fopen(filename,'r');
  % if standard reading is not the correct format for a given binary
  % file, activate the following:
  %fid = fopen(filename,'r','ieee-be');
  
  magicstr = char(fread(fid,3,'char'))';
  if ~isequal(magicstr,'DEV')
    error('read magicstr doe not indicate Dune Sparse Matrix!');
  end;

  magicint = fread(fid,1,'int');
  magicdouble = fread(fid,1,'double');
  
  if (magicint~=111) | (magicdouble~=111.0)
    error(['magic numbers not read correctly, change the binary format in' ...
	   ' this reading routine!!']);
  end;
  
  ncols = fread(fid,1,'int');
  nrows = fread(fid,1,'int');
  
  disp(['generating ',num2str(nrows),'x',num2str(ncols),' matrix.']);
  v = fread(fid,nrows*ncols,'double');
  A = zeros(nrows,ncols);
  
  for i=1:nrows
      A(i,:) = v((i-1)*ncols + 1 : i*ncols);
  end; 
  
  eofstr = char(fread(fid,3,'char'))';
  if ~isequal(eofstr,'EOF')
    error('read eofstr does not indicate end of binary file!');
  end;

  fclose(fid);
