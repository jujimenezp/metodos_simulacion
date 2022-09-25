set g
set xl 'x'
set yl 'y'
set size 1,1

set term pdf
set o 'data/Taller1/punto_3b.pdf'

set title 'Jupiter'
plot 'data/Taller1/punto_3b.dat' u 1:2 w lp lc 'blue' notitle

set title 'Sol'
plot 'data/Taller1/punto_3b.dat' u 3:4 w l lc 'red' notitle
