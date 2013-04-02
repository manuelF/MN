#!/usr/local/bin/octave -qf
source("leer_datos.m");
source("dibujarHistyAjuste.m");
source("dibujarHist.m");
source("GGDpdf_c.m");
args = argv();
d = leer_datos(args{1});

figure("visible", "off");
hold on
disp(args);

dibujarHist(d);
dibujarAjuste(d, str2double(args{2}),
                      str2double(args{3}),
                      str2double(args{4}),
                      0.0,0.0,0.0);

if(size(args,1)==8)
dibujarAjuste(d, str2double(args{5}),
                      str2double(args{6}),
                      str2double(args{7}),
                      1.0,0.0,0.0);


	fname =args{8};
	print(fname, '-dpng');
else
	if(size(args,1)==11)
	dibujarAjuste(d, str2double(args{5}),
	                  str2double(args{6}),
	                  str2double(args{7}),
	                  1.0,0.0,0.0);
	dibujarAjuste(d, str2double(args{8}),
	                 str2double(args{9}),
	                 str2double(args{10}),
	                 0.0,0.0,1.0);


		fname =args{11};
		print(fname, '-dpng');
	else


		if(size(args,1)==5)
			fname =args{5};
			print(fname, '-dpng');
		else
			print("plot.png", '-dpng');
		end
	end

end