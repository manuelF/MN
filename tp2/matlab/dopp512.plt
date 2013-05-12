set terminal pngcairo enhanced font 'Verdana,10'
set output 'dopp512.png'
set autoscale
set title 'Dopp 512'
set xlabel 'Punto'
set ylabel 'Valor'
set key top right

set style line 1 lt 1 lw 1 pt 1 linecolor rgb "red"
set style line 2 lt 1 lw 1 pt 1 linecolor rgb "blue"
set style line 3 lt 1 lw 1 pt 1 linecolor rgb "green"
plot		"orig1"  using 1 title 'Original'  with lines ls 1, \
	"mod1"  using 1 title 'Modified'  with lines ls 2, \
	"rec1"  using 1 title 'Recovered'  with lines ls 3

