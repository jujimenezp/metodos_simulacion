set g
set xl 'x'
set yl 'y'
set size 1,1

set term pdf
set o 'data/Taller1/punto_3d.pdf'

set title 'Jupiter'
plot 'data/Taller1/punto_3d.dat' u 1:2 w lp lc 'blue' notitle

set title 'Sol'
plot 'data/Taller1/punto_3d.dat' u 3:4 w l lc 'black' notitle

set title 'Troyano'
plot 'data/Taller1/punto_3d.dat' u 5:6 w l lc 'red' notitle
