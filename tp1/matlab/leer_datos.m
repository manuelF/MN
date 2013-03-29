function [data, n] = leer_datos(fname)
%
%	[data, n] = leer_datos(fname)
%	Lee datos desde archivo fname. 
%     Formato:
%     L1: n
%     L2: X1 X2 X3 ... Xn
%

fp = fopen(fname);
if fp == -1
	disp(['Error: no se puede abrir archivo ' fname]);
	return;
end	

% lee cantidad de datos
n = fscanf(fp,'%d\n',1);

% datos
data = fscanf(fp, '%f');

fclose(fp);


