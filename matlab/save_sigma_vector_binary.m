function save_sigma_vector_binary(filename, sigma, id)

  fid = fopen(filename,'w');
  % if standard reading is not the correct format for a given binary
  % file, activate the following:
  %fid = fopen(filename,'r','ieee-be');
  
  fwrite(fid,'DDV','char*1');
  
  fwrite(fid,111,'int');
  fwrite(fid,111.0,'double');
  
  fwrite(fid,length(sigma),'int');

  zippervector = zeros(1,2*length(sigma));
  zippervector(1:2:end) = id;
  zippervector(2:2:end) = sigma;
  fwrite(fid,zippervector,'double');

  fwrite(fid,'EOF','char*1');
  
  fclose(fid);