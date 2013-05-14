set terminal pngcairo enhanced font 'Verdana,10' size 1024,768
set output 'ramp1234.png'
set autoscale
set title 'Ramp 1234'
set xlabel 'Punto'
set ylabel 'Valor'
set key top right

set style line 1 lt 1 lw 1 pt 1 linecolor rgb "red"
set style line 2 lt 1 lw 1 pt 1 linecolor rgb "blue"
set style line 3 lt 1 lw 1 pt 1 linecolor rgb "green"
plot		"orig2"  using 1 title 'Original'  with lines ls 1, \
	"mod2"  using 1 title 'Modified'  with lines ls 2, \
	"rec2"  using 1 title 'Recovered'  with lines ls 3

