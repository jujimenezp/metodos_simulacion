set grid
set term jpeg
set o 'data/1b.jpg'

set xl 'Tiempo'
set yl 'Posición x'
plot 'data/1b.dat' w l notitle
