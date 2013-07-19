
set output 'g450.png'

set terminal pngcairo enhanced font 'Verdana,10' size 1024,768
set autoscale
##set xlabel 'Numero de muestra'
set ylabel 'Magnitud'

set style line 1 lt 1 lw 1 pt 1 linecolor rgb "red"
set style line 2 lt 1 lw 1 pt 1 linecolor rgb "blue"
set style line 3 lt 1 lw 1 pt 1 linecolor rgb "green"


set tmargin 0
set bmargin 0
set lmargin 10
set rmargin 3
unset xtics

##set multiplot layout 4,1 title "G450 con ruido gaussiano (varianza 10) y filtro exponencial"
set multiplot layout 4,1 title "G450 con ruido salt and pepper y filtro exponencial"

set key autotitle column nobox samplen 1 noenhanced
unset title
set style data boxes

plot "orig4"  using 1 title 'Original'  with lines ls 1
plot "mod4"  using 1 title 'Modificada'  with lines ls 2
##set xtics nomirror
##set tics scale 0 font ",8"
set xlabel "Descomposicion en el espectro de frecuencias"
plot "rec4"  using 1 title 'Recuperada'  with lines ls 3

unset multiplot
