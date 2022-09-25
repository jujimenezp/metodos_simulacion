set g
set xl 't'
set yl 'x(t)'


set term pdf
set o 'data/Taller1/punto_3e.pdf'

set title 'Troyano'
do for [i=1:15]{
set arrow from 6157.46*i+1000,460 to 6157.46*i+1000,540 nohead
}
do for [i=1:4]{
set arrow from 76440,460 to 76440,540 nohead lc rgb 'red'
}

plot 'data/Taller1/punto_3e.dat' u 1:2 w l lc 'red' notitle
