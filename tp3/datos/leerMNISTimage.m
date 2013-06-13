function img = leerMNISTimage(fname, limit)
% All the integers in the files are stored in the MSB first (high/big endian) 
% TEST SET IMAGE FILE (t10k-images-idx3-ubyte):
% 
% [offset] [type]          [value]          [description] 
% 0000     32 bit integer  0x00000803(2051) magic number 
% 0004     32 bit integer  60000            number of images 
% 0008     32 bit integer  28               number of rows 
% 0012     32 bit integer  28               number of columns 
% 0016     unsigned byte   ??               pixel 
% 0017     unsigned byte   ??               pixel 
% ........ 
% xxxx     unsigned byte   ??               pixel
% Pixels are organized row-wise. Pixel values are 0 to 255. 0 means background, 255 means foreground. 


  fp = fopen(fname);
  if fp == -1
   error(['Error: no se puede abrir archivo ' fname]);
  end	
  
  magicnumber = fread(fp,1,'uint32','ieee-be');
  nimgs = fread(fp,1,'uint32','ieee-be');
  nrows = fread(fp,1,'uint32','ieee-be');
  ncols = fread(fp,1,'uint32','ieee-be');
  toread=min(nimgs,limit);
  img = zeros(nrows*ncols,toread,'uint8');
  
    
  for imind=1:toread,
     im = fread(fp,[ncols*nrows],'*uint8');
     img(:,imind) = im';
  end
  
  fclose(fp);
