set grid

set xl 'Tiempo'
set yl 'Casos'

set term pdf

set o 'data/Taller1/punto_1a.pdf'
plot 'data/Taller1/punto_1a.dat' u 1:2 w l lc rgb 'blue' title 's(t)',\
     '' u 1:3 w l lc rgb 'red' title 'i(t)',\
     '' u 1:4 w l lc rgb 'green' title 'r(t)'
