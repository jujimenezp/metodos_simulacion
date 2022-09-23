set ytics 0.1
set xtics 1
set grid

set xl 'r'
set yl 'R(r)'

set term pdf
set o 'data/Taller1/punto_2c.pdf'

plot 'data/Taller1/punto_2c.dat' u 1:2 w lp lc rgb 'blue' pt 7 ps 0.3 title 'Calculado',\
     'data/Taller1/punto_2d.dat' u 1:2 w l lc rgb 'black' title 'Teórico'
