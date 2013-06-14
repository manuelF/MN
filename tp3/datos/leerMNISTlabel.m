function labels = leerMNISTlabel(fname, from, limit)
% All the integers in the files are stored in the MSB first (high/big endian) 
% [offset] [type]          [value]          [description] 
% 0000     32 bit integer  0x00000801(2049) magic number (MSB first) 
% 0004     32 bit integer  60000            number of items 
% 0008     unsigned byte   ??               label 
% 0009     unsigned byte   ??               label 
% ........ 
% xxxx     unsigned byte   ??               label
% The labels values are 0 to 9.


  fp = fopen(fname,'r');
  if fp == -1
   error(['Error: no se puede abrir archivo ' fname]);
  end	
  
  magicnumber = fread(fp,1,'uint32','ieee-be');
  nitems = fread(fp,1,'uint32','ieee-be');
  nitems = min(nitems,limit);
  
  labels2 = fread(fp,nitems,'*uint8');
  labels = labels2(from:limit);
  fclose(fp);
