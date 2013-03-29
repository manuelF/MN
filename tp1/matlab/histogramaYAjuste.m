#!/usr/local/bin/octave -qf
source("leer_datos.m")
source("dibujarHistyAjuste.m")
source("GGDpdf_c.m")
args = argv();
d = leer_datos(args{1});
dibujarHistyAjuste(d, str2double(args{2}),
                      str2double(args{3}),
                      str2double(args{4}))
print -dpng plot.png
