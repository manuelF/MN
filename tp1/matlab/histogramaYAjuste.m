#!/usr/local/bin/octave -qf
source("leer_datos.m");
source("dibujarHistyAjuste.m");
source("dibujarHist.m");
source("GGDpdf_c.m");
args = argv();
d = leer_datos(args{1});

figure;
hold on
disp(args);

dibujarHist(d);
dibujarAjuste(d, str2double(args{2}),
                      str2double(args{3}),
                      str2double(args{4}),
                      'b');
if(size(args,1)==8)
dibujarAjuste(d, str2double(args{5}),
                      str2double(args{6}),
                      str2double(args{7}),
                      'r');


	fname =args{8};
	print(fname, '-dpng');
else

	#print -dpng plot.png

end