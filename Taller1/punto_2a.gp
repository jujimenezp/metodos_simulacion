set ytics 0.1
set xtics 0.5
set grid

set xl 'r'
set yl 'R(r)'

set term pdf
set o 'data/Taller1/punto_2a.pdf'

plot 'data/Taller1/punto_2a.dat' u 1:2 w l lc rgb 'blue' notitle
