set grid
set term jpeg
set o 'data/1_2000.jpg'

set xl 'x (celdas)'
set yl 'rho'

plot 'data/1_2000.dat' u 1:3  w l title '2000 pasos'

set o 'data/1_6000.jpg'
plot 'data/1_6000.dat' u 1:3  w l title '6000 pasos'
