
set output 'dopp512.png'
set title 'Dopp 512'

set terminal pngcairo enhanced font 'Verdana,10' size 1024,768
set autoscale
set xlabel 'Numero de muestra'
set ylabel 'Magnitud'
set key top right

set style line 1 lt 1 lw 1 pt 1 linecolor rgb "red"
set style line 2 lt 1 lw 1 pt 1 linecolor rgb "blue"
set style line 3 lt 1 lw 1 pt 1 linecolor rgb "green"
plot		"orig1"  using 1 title 'Original'  with lines ls 1, \
	"mod1"  using 1 title 'Modificada'  with lines ls 2, \
	"rec1"  using 1 title 'Recuperada'  with lines ls 3

