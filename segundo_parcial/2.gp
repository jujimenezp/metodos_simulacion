set grid
set term jpeg
set o 'data/2_envolventes.jpg'

set xl 'x (celdas)'
set yl 'rho'
set yrange [-1:1]
set ytics 0.1

plot 'data/2_envolventes.dat' u 1:2 w l title 'rho mínimo', 'data/2_envolventes.dat' u 1:3 w l title 'rho maximo'

set o "data/2_ondas.jpeg"

set ytics 0.1
plot for [i=80000:82180:20] 'data/2_'.i.'.dat' u 1:3 w l notitle, 'data/2_envolventes.dat' u 1:2 title 'rho mínimo', 'data/2_envolventes.dat' u 1:3 title 'rho maximo'
