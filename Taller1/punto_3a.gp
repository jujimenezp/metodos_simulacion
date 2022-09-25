set g
set xl 'x'
set yl 'y'
set size 1,1

set term pdf
set o 'data/Taller1/punto_3a.pdf'

set title 'Jupiter'
plot 'data/Taller1/punto_3a.dat' u 1:2 w l lc 'blue' notitle

set title 'Sol'
plot 'data/Taller1/punto_3a.dat' u 3:4 w l lc 'red' notitle
